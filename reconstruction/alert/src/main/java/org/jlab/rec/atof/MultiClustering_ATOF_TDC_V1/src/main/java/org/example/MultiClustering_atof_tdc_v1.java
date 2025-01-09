
package org.jlab.rec.atof.MultiClustering_ATOF_TDC_V1;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.io.HipoWriter;

import java.util.ArrayList;
import java.util.List;

public class MultiClustering_atof_tdc_v1 {

    private static final double Z_THRESHOLD = 200.0;
    private static final double PHI_THRESHOLD = 0.2;
    private static final double TIME_THRESHOLD = 2.5;
    private static final double VEFF = 200.0;
    private static final int NUM_BARS = 60;
    private static final int WEDGES_PER_BAR = 10;
    private static final double WEDGE_SPACING = 30.0;
    private static final double TDC_TO_NS = 0.015625; // Conversion factor for TDC to nanoseconds

    public static void main(String[] args) {
        if (args.length < 3) {
            System.err.println("Usage: java MultiClustering <input.hipo> <output.hipo> <schema.json>");
            System.exit(1);
        }

        String inputHipoFile = args[0];
        String outputHipoFile = args[1];
        String schemaJsonFile = args[2];

        HipoReader reader = new HipoReader();
        reader.open(inputHipoFile);

        SchemaFactory schemaFactory = reader.getSchemaFactory();
        try {
            schemaFactory.readFile(schemaJsonFile);
        } catch (Exception e) {
            System.err.println("Error loading schema: " + e.getMessage());
            System.exit(1);
        }

        Schema recSchema = schemaFactory.getSchema("ATOF::rec");
        Bank tdcBank = new Bank(schemaFactory.getSchema("ATOF::tdc"));
        Event event = new Event();

        HipoWriter writer = new HipoWriter();
        writer.getSchemaFactory().addSchema(recSchema);
        writer.open(outputHipoFile);

        int eventId = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(tdcBank);

            System.out.printf("\nProcessing Event ID: %d with %d hits...\n", eventId, tdcBank.getRows());

            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();
            extractHits(tdcBank, barHits, wedgeHits);

            System.out.println("Hits in this event:");
            printHits(barHits, "Bar");
            printHits(wedgeHits, "Wedge");

            List<Cluster> barClusters = formBarClusters(barHits);
            List<Cluster> barWedgeClusters = formBarWedgeClusters(barClusters, wedgeHits);

            Bank recBank = new Bank(recSchema, barClusters.size() + barWedgeClusters.size());
            int clusterId = 0;

            System.out.println("Bar-only Clusters:");
            for (Cluster cluster : barClusters) {
                storeCluster(recBank, clusterId++, cluster);
                printCluster(cluster, eventId, "Bar Only");
            }

            System.out.println("Bar + Wedge Clusters:");
            for (Cluster cluster : barWedgeClusters) {
                storeCluster(recBank, clusterId++, cluster);
                printCluster(cluster, eventId, "Bar + Wedge");
            }

            event.write(recBank);
            writer.addEvent(event);
            eventId++;
        }

        reader.close();
        writer.close();
    }

    private static void extractHits(Bank tdcBank, List<Hit> barHits, List<Hit> wedgeHits) {
        for (int i = 0; i < tdcBank.getRows(); i++) {
            int sector = tdcBank.getByte("sector", i);
            int layer = tdcBank.getByte("layer", i);
            int component = tdcBank.getShort("component", i);
            int order = tdcBank.getByte("order", i);
            int tdc = tdcBank.getInt("TDC", i);
            int tot = tdcBank.getInt("ToT", i);
            double time = tdc * TDC_TO_NS;
            double phi = -Math.PI + (2 * Math.PI * component) / NUM_BARS;
            double z = (layer == 0) ? 0.0 : ((component % WEDGES_PER_BAR) - 4.5) * WEDGE_SPACING;

            Hit hit = new Hit(sector, layer, component, order, tdc, tot, time, z, phi);

            if (layer == 0) barHits.add(hit);
            else wedgeHits.add(hit);
        }
    }

    private static List<Cluster> formBarClusters(List<Hit> barHits) {
        List<Cluster> clusters = new ArrayList<>();
        for (int i = 0; i < barHits.size(); i++) {
            for (int j = i + 1; j < barHits.size(); j++) {
                Hit hit1 = barHits.get(i);
                Hit hit2 = barHits.get(j);
                if (hit1.component == hit2.component) {
                    double zBar = (VEFF / 2.0) * (hit2.time - hit1.time);
                    double tMin = Math.min(hit1.time, hit2.time);
                    double energy = hit1.tot + hit2.tot;
                    Cluster cluster = new Cluster(zBar, hit1.phi, tMin, energy);
                    cluster.hits.add(hit1);
                    cluster.hits.add(hit2);
                    clusters.add(cluster);
                }
            }
        }
        return clusters;
    }

    private static List<Cluster> formBarWedgeClusters(List<Cluster> barClusters, List<Hit> wedgeHits) {
        List<Cluster> clusters = new ArrayList<>();
        for (Cluster barCluster : barClusters) {
            for (Hit wedgeHit : wedgeHits) {
                double deltaZ = Math.abs(barCluster.z - wedgeHit.z);
                double deltaPhi = Math.abs(barCluster.phi - wedgeHit.phi);
                double deltaTime = Math.abs(barCluster.time - wedgeHit.time);

                if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                    double energy = barCluster.energy + wedgeHit.tot;
                    Cluster cluster = new Cluster(barCluster.z, barCluster.phi, Math.min(barCluster.time, wedgeHit.time), energy);
                    cluster.hits.addAll(barCluster.hits);
                    cluster.hits.add(wedgeHit);
                    cluster.deltaZ = deltaZ;
                    cluster.deltaPhi = deltaPhi;
                    cluster.deltaTime = deltaTime;
                    clusters.add(cluster);
                }
            }
        }
        return clusters;
    }

    private static void printHits(List<Hit> hits, String type) {
        System.out.printf("  %s Hits (%d total):\n", type, hits.size());
        for (Hit hit : hits) {
            System.out.printf("    Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, TDC: %d, ToT: %d, Time: %.2f ns, Z: %.2f mm, Phi: %.2f rad\n",
                    hit.sector, hit.layer, hit.component, hit.order, hit.tdc, hit.tot, hit.time, hit.z, hit.phi);
        }
    }

    private static void printCluster(Cluster cluster, int eventId, String type) {
        System.out.printf("%s Cluster ID: %d, Event ID: %d -> Z: %.2f mm, Phi: %.2f rad, Time: %.2f ns, Energy: %.2f MeV, Size: %d\n",
                type, cluster.id, eventId, cluster.z, cluster.phi, cluster.time, cluster.energy, cluster.hits.size());
        for (Hit hit : cluster.hits) {
            System.out.printf("      Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, TDC: %d, ToT: %d, Time: %.2f ns, Z: %.2f mm, Phi: %.2f rad\n",
                    hit.sector, hit.layer, hit.component, hit.order, hit.tdc, hit.tot, hit.time, hit.z, hit.phi);
        }
        if ("Bar + Wedge".equals(type)) {
            System.out.printf("    Thresholds -> Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n",
                    cluster.deltaZ, cluster.deltaPhi, cluster.deltaTime);
        }
    }

    private static void storeCluster(Bank recBank, int clusterId, Cluster cluster) {
        recBank.putShort("id", clusterId, (short) clusterId);
        recBank.putShort("nhits", clusterId, (short) cluster.hits.size());
        recBank.putFloat("z", clusterId, (float) cluster.z);
        recBank.putFloat("phi", clusterId, (float) cluster.phi);
        recBank.putFloat("time", clusterId, (float) cluster.time);
        recBank.putFloat("energy", clusterId, (float) cluster.energy);
    }

    static class Hit {
        int sector, layer, component, order, tdc, tot;
        double time, z, phi;

        Hit(int sector, int layer, int component, int order, int tdc, int tot, double time, double z, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.tdc = tdc;
            this.tot = tot;
            this.time = time;
            this.z = z;
            this.phi = phi;
        }
    }

    static class Cluster {
        double z, phi, time, deltaZ, deltaPhi, deltaTime, energy;
        List<Hit> hits = new ArrayList<>();
        int id;

        Cluster(double z, double phi, double time, double energy) {
            this.z = z;
            this.phi = phi;
            this.time = time;
            this.energy = energy;
        }
    }
}






/*
//USE OF TDC BANK

package org.jlab.rec.atof.MultiClustering_ATOF_TDC_V1;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.io.HipoWriter;

import java.util.ArrayList;
import java.util.List;

public class MultiClustering_atof_tdc_v1 {

    private static final double Z_THRESHOLD = 200.0;
    private static final double PHI_THRESHOLD = 0.2;
    private static final double TIME_THRESHOLD = 2.5;
    private static final double VEFF = 200.0;
    private static final int NUM_BARS = 60;
    private static final int WEDGES_PER_BAR = 10;
    private static final double WEDGE_SPACING = 30.0;
    private static final double TDC_TO_NS = 0.015625; // Conversion factor for TDC to nanoseconds

    public static void main(String[] args) {
        if (args.length < 3) {
            System.err.println("Usage: java MultiClustering <input.hipo> <output.hipo> <schema.json>");
            System.exit(1);
        }

        String inputHipoFile = args[0];
        String outputHipoFile = args[1];
        String schemaJsonFile = args[2];

        HipoReader reader = new HipoReader();
        reader.open(inputHipoFile);

        SchemaFactory schemaFactory = reader.getSchemaFactory();
        try {
            schemaFactory.readFile(schemaJsonFile);
        } catch (Exception e) {
            System.err.println("Error loading schema: " + e.getMessage());
            System.exit(1);
        }

        Schema recSchema = schemaFactory.getSchema("ATOF::rec");
        Bank tdcBank = new Bank(schemaFactory.getSchema("ATOF::tdc"));
        Event event = new Event();

        HipoWriter writer = new HipoWriter();
        writer.getSchemaFactory().addSchema(recSchema);
        writer.open(outputHipoFile);

        int eventId = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(tdcBank);

            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();
            extractHits(tdcBank, barHits, wedgeHits);

            List<Cluster> barClusters = formBarClusters(barHits);
            List<Cluster> barWedgeClusters = formBarWedgeClusters(barClusters, wedgeHits);

            Bank recBank = new Bank(recSchema, barClusters.size() + barWedgeClusters.size());
            int clusterId = 0;

            System.out.printf("Processing Event ID: %d\n", eventId);

            for (Cluster cluster : barClusters) {
                storeCluster(recBank, clusterId++, cluster);
            }

            for (Cluster cluster : barWedgeClusters) {
                storeCluster(recBank, clusterId++, cluster);
            }

            event.write(recBank);
            writer.addEvent(event);
            eventId++;
        }

        reader.close();
        writer.close();
    }

    private static void extractHits(Bank tdcBank, List<Hit> barHits, List<Hit> wedgeHits) {
        for (int i = 0; i < tdcBank.getRows(); i++) {
            int sector = tdcBank.getByte("sector", i);
            int layer = tdcBank.getByte("layer", i);
            int component = tdcBank.getShort("component", i);
            int order = tdcBank.getByte("order", i);
            int tdc = tdcBank.getInt("TDC", i);
            int tot = tdcBank.getInt("ToT", i);
            double time = tdc * TDC_TO_NS;
            double phi = -Math.PI + (2 * Math.PI * component) / NUM_BARS;
            double z = (layer == 0) ? 0.0 : ((component % WEDGES_PER_BAR) - 4.5) * WEDGE_SPACING;

            Hit hit = new Hit(sector, layer, component, order, tdc, tot, time, z, phi);

            if (layer == 0) barHits.add(hit);
            else wedgeHits.add(hit);
        }
    }

    private static List<Cluster> formBarClusters(List<Hit> barHits) {
        List<Cluster> clusters = new ArrayList<>();
        for (int i = 0; i < barHits.size(); i++) {
            for (int j = i + 1; j < barHits.size(); j++) {
                Hit hit1 = barHits.get(i);
                Hit hit2 = barHits.get(j);
                if (hit1.component == hit2.component) {
                    double zBar = (VEFF / 2.0) * (hit2.time - hit1.time);
                    double tMin = Math.min(hit1.time, hit2.time);
                    double energy = hit1.tot + hit2.tot;
                    Cluster cluster = new Cluster(zBar, hit1.phi, tMin, energy);
                    cluster.hits.add(hit1);
                    cluster.hits.add(hit2);
                    clusters.add(cluster);
                }
            }
        }
        return clusters;
    }

    private static List<Cluster> formBarWedgeClusters(List<Cluster> barClusters, List<Hit> wedgeHits) {
        List<Cluster> clusters = new ArrayList<>();
        for (Cluster barCluster : barClusters) {
            for (Hit wedgeHit : wedgeHits) {
                double deltaZ = Math.abs(barCluster.z - wedgeHit.z);
                double deltaPhi = Math.abs(barCluster.phi - wedgeHit.phi);
                double deltaTime = Math.abs(barCluster.time - wedgeHit.time);

                if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                    double energy = barCluster.energy + wedgeHit.tot;
                    Cluster cluster = new Cluster(barCluster.z, barCluster.phi, Math.min(barCluster.time, wedgeHit.time), energy);
                    cluster.hits.addAll(barCluster.hits);
                    cluster.hits.add(wedgeHit);
                    cluster.deltaZ = deltaZ;
                    cluster.deltaPhi = deltaPhi;
                    cluster.deltaTime = deltaTime;
                    clusters.add(cluster);
                }
            }
        }
        return clusters;
    }

    private static void storeCluster(Bank recBank, int clusterId, Cluster cluster) {
        recBank.putShort("id", clusterId, (short) clusterId);
        recBank.putShort("nhits", clusterId, (short) cluster.hits.size());
        recBank.putFloat("z", clusterId, (float) cluster.z);
        recBank.putFloat("phi", clusterId, (float) cluster.phi);
        recBank.putFloat("time", clusterId, (float) cluster.time);
        recBank.putFloat("energy", clusterId, (float) cluster.energy);
    }

    static class Hit {
        int sector, layer, component, order, tdc, tot;
        double time, z, phi;

        Hit(int sector, int layer, int component, int order, int tdc, int tot, double time, double z, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.tdc = tdc;
            this.tot = tot;
            this.time = time;
            this.z = z;
            this.phi = phi;
        }
    }

    static class Cluster {
        double z, phi, time, deltaZ, deltaPhi, deltaTime, energy;
        List<Hit> hits = new ArrayList<>();

        Cluster(double z, double phi, double time, double energy) {
            this.z = z;
            this.phi = phi;
            this.time = time;
            this.energy = energy;
        }
    }
}
*/












/*
//This performs multiclustering, pushes cluster output to ATOF::rec bank, bank is created and output is pushed there!!

package org.jlab.rec.atof.MultiClustering_ATOF_V1;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.io.HipoWriter;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class MultiClustering_atof_v1 {

    private static final double Z_THRESHOLD = 200.0;
    private static final double PHI_THRESHOLD = 0.2;
    private static final double TIME_THRESHOLD = 2.5;
    private static final double VEFF = 200.0;
    private static final int NUM_BARS = 60;
    private static final int WEDGES_PER_BAR = 10;
    private static final double WEDGE_SPACING = 30.0;

    public static void main(String[] args) {
        if (args.length < 3) {
            System.err.println("Usage: java MultiClustering <input.hipo> <output.hipo> <schema.json>");
            System.exit(1);
        }

        String inputHipoFile = args[0];
        String outputHipoFile = args[1];
        String schemaJsonFile = args[2];

        HipoReader reader = new HipoReader();
        reader.open(inputHipoFile);

        SchemaFactory schemaFactory = reader.getSchemaFactory();
        try {
            schemaFactory.readFile(schemaJsonFile);
        } catch (Exception e) {
            System.err.println("Error loading schema: " + e.getMessage());
            System.exit(1);
        }

        Schema recSchema = schemaFactory.getSchema("ATOF::rec");
        Bank adcBank = new Bank(schemaFactory.getSchema("ATOF::adc"));
        Event event = new Event();

        HipoWriter writer = new HipoWriter();
        writer.getSchemaFactory().addSchema(recSchema);
        writer.open(outputHipoFile);

        int eventId = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(adcBank);

            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();
            extractHits(adcBank, barHits, wedgeHits);

            List<Cluster> barClusters = formBarClusters(barHits);
            List<Cluster> barWedgeClusters = formBarWedgeClusters(barClusters, wedgeHits);

            Bank recBank = new Bank(recSchema, barClusters.size() + barWedgeClusters.size());
            int clusterId = 0;

            System.out.printf("Processing Event ID: %d\n", eventId);

            System.out.println("Bar-only Clusters:");
            for (Cluster cluster : barClusters) {
                storeCluster(recBank, clusterId++, cluster);
                printCluster(cluster, eventId, "Bar Only");
            }

            System.out.println("Bar + Wedge Clusters:");
            for (Cluster cluster : barWedgeClusters) {
                storeCluster(recBank, clusterId++, cluster);
                printCluster(cluster, eventId, "Bar + Wedge");
            }

            event.write(recBank);
            writer.addEvent(event);
            eventId++;
        }

        reader.close();
        writer.close();
    }

    private static void extractHits(Bank adcBank, List<Hit> barHits, List<Hit> wedgeHits) {
        for (int i = 0; i < adcBank.getRows(); i++) {
            int sector = adcBank.getByte("sector", i);
            int layer = adcBank.getByte("layer", i);
            int component = adcBank.getShort("component", i);
            int order = adcBank.getByte("order", i);
            int adc = adcBank.getInt("ADC", i);
            float time = adcBank.getFloat("time", i);
            int ped = adcBank.getShort("ped", i);
            double phi = -Math.PI + (2 * Math.PI * component) / NUM_BARS;
            double z = (layer == 0) ? 0.0 : ((component % WEDGES_PER_BAR) - 4.5) * WEDGE_SPACING;
            Hit hit = new Hit(sector, layer, component, order, adc, time, ped, z, phi);

            if (layer == 0) barHits.add(hit);
            else wedgeHits.add(hit);
        }
    }

    private static List<Cluster> formBarClusters(List<Hit> barHits) {
        List<Cluster> clusters = new ArrayList<>();
        for (int i = 0; i < barHits.size(); i++) {
            for (int j = i + 1; j < barHits.size(); j++) {
                Hit hit1 = barHits.get(i);
                Hit hit2 = barHits.get(j);
                if (hit1.component == hit2.component) {
                    double zBar = (VEFF / 2.0) * (hit2.time - hit1.time);
                    double tMin = Math.min(hit1.time, hit2.time);
                    double energy = hit1.adc + hit2.adc;
                    Cluster cluster = new Cluster(zBar, hit1.phi, tMin, energy);
                    cluster.hits.add(hit1);
                    cluster.hits.add(hit2);
                    clusters.add(cluster);
                }
            }
        }
        return clusters;
    }

    private static List<Cluster> formBarWedgeClusters(List<Cluster> barClusters, List<Hit> wedgeHits) {
        List<Cluster> clusters = new ArrayList<>();
        for (Cluster barCluster : barClusters) {
            for (Hit wedgeHit : wedgeHits) {
                double deltaZ = Math.abs(barCluster.z - wedgeHit.z);
                double deltaPhi = Math.abs(barCluster.phi - wedgeHit.phi);
                double deltaTime = Math.abs(barCluster.time - wedgeHit.time);

                if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                    double energy = barCluster.energy + wedgeHit.adc;
                    Cluster cluster = new Cluster(barCluster.z, barCluster.phi, Math.min(barCluster.time, wedgeHit.time), energy);
                    cluster.hits.addAll(barCluster.hits);
                    cluster.hits.add(wedgeHit);
                    cluster.deltaZ = deltaZ;
                    cluster.deltaPhi = deltaPhi;
                    cluster.deltaTime = deltaTime;
                    clusters.add(cluster);
                }
            }
        }
        return clusters;
    }

    private static void storeCluster(Bank recBank, int clusterId, Cluster cluster) {
        recBank.putShort("id", clusterId, (short) clusterId);
        recBank.putShort("nhits", clusterId, (short) cluster.hits.size());
        recBank.putFloat("z", clusterId, (float) cluster.z);
        recBank.putFloat("phi", clusterId, (float) cluster.phi);
        recBank.putFloat("time", clusterId, (float) cluster.time);
        recBank.putFloat("energy", clusterId, (float) cluster.energy);
    }

    private static void printCluster(Cluster cluster, int eventId, String type) {
        System.out.printf("%s Cluster ID: %d, Event ID: %d -> Z: %.2f mm, Phi: %.2f rad, Time: %.2f ns, Energy: %.2f MeV, Size: %d\n",
                type, cluster.id, eventId, cluster.z, cluster.phi, cluster.time, cluster.energy, cluster.hits.size());
        for (Hit hit : cluster.hits) {
            if (type.equals("Bar Only") || (hit.layer == 0 && type.equals("Bar + Wedge"))) {
                System.out.printf("  Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                        hit.sector, hit.layer, hit.component, hit.order, hit.adc, hit.time, hit.ped, hit.phi);
            } else {
                System.out.printf("  Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Z: %.2f mm, Phi: %.2f rad\n",
                        hit.sector, hit.layer, hit.component, hit.order, hit.adc, hit.time, hit.ped, hit.z, hit.phi);
            }
        }
        if ("Bar + Wedge".equals(type)) {
            System.out.printf("    Thresholds -> Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n",
                    cluster.deltaZ, cluster.deltaPhi, cluster.deltaTime);
        }
    }

    static class Hit {
        int sector, layer, component, order, adc, ped;
        double time, z, phi;

        Hit(int sector, int layer, int component, int order, int adc, double time, int ped, double z, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.time = time;
            this.ped = ped;
            this.z = z;
            this.phi = phi;
        }
    }

    static class Cluster {
        double z, phi, time, deltaZ, deltaPhi, deltaTime, energy;
        List<Hit> hits = new ArrayList<>();
        int id;

        Cluster(double z, double phi, double time, double energy) {
            this.z = z;
            this.phi = phi;
            this.time = time;
            this.energy = energy;
        }
    }
}

*/
