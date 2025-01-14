//This is the code for bar and wedge multiclustering. Per event multiple cluster exists. Code prints the cluster output. 

package org.jlab.rec.atof.MultiClustering_ATOF_TDC_V2;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.io.HipoWriter;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class MultiClustering_atof_tdc_v2 {

    private static final double Z_THRESHOLD = 300.0;
    private static final double PHI_THRESHOLD = 0.3;
    private static final double TIME_THRESHOLD = 3.0;
    private static final double VEFF = 200.0;
    private static final int NUM_BARS = 60;
    private static final int WEDGES_PER_BAR = 10;
    private static final double WEDGE_SPACING = 30.0;
    private static final double TDC_TO_NS = 0.015625;

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

            System.out.println("\nHits in this event:");
            printHits(barHits, "Bar");
            printHits(wedgeHits, "Wedge");

            List<Cluster> barClusters = formBarClusters(barHits);
            List<Cluster> barWedgeClusters = formBarWedgeClusters(barClusters, wedgeHits);

            System.out.println("\nBar-only Clusters:");
            printClusters(barClusters, "Bar-Only");

            System.out.println("\nBar+Wedge Clusters:");
            printClusters(barWedgeClusters, "Bar+Wedge");

            Bank recBank = new Bank(recSchema, barClusters.size() + barWedgeClusters.size());
            storeClusters(recBank, barClusters);
            storeClusters(recBank, barWedgeClusters);

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
        Map<Integer, List<Hit>> hitsByComponent = new HashMap<>();

        for (Hit hit : barHits) {
            hitsByComponent.computeIfAbsent(hit.component, k -> new ArrayList<>()).add(hit);
        }

        for (Map.Entry<Integer, List<Hit>> entry : hitsByComponent.entrySet()) {
            List<Hit> hits = entry.getValue();

            Hit leftHit = null;
            Hit rightHit = null;

            for (Hit hit : hits) {
                if (hit.order == 0) leftHit = hit;
                if (hit.order == 1) rightHit = hit;
            }

            if (leftHit != null && rightHit != null) {
                double z = (VEFF / 2.0) * (rightHit.time - leftHit.time);
                double tMin = Math.min(leftHit.time, rightHit.time);
                double energy = leftHit.tot + rightHit.tot;

                Cluster cluster = new Cluster(z, leftHit.phi, tMin, energy, false); // Exclude Z for Bar clusters
                cluster.hits.add(leftHit);
                cluster.hits.add(rightHit);

                clusters.add(cluster);
            }
        }

        return clusters;
    }

    private static List<Cluster> formBarWedgeClusters(List<Cluster> barClusters, List<Hit> wedgeHits) {
        List<Cluster> clusters = new ArrayList<>();

        for (Cluster barCluster : barClusters) {
            Cluster combinedCluster = new Cluster(barCluster.z, barCluster.phi, barCluster.time, barCluster.energy, true); // Include Z for Bar+Wedge clusters
            combinedCluster.hits.addAll(barCluster.hits);

            for (Hit wedgeHit : wedgeHits) {
                if (wedgeHit.component / WEDGES_PER_BAR == barCluster.hits.get(0).component / WEDGES_PER_BAR) {
                    double deltaZ = Math.abs(barCluster.z - wedgeHit.z);
                    double deltaPhi = Math.abs(barCluster.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(barCluster.time - wedgeHit.time);

                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        combinedCluster.hits.add(wedgeHit);
                        combinedCluster.energy += wedgeHit.tot;
                        combinedCluster.time = Math.min(combinedCluster.time, wedgeHit.time);
                    }
                }
            }

            if (combinedCluster.hits.size() > barCluster.hits.size()) {
                clusters.add(combinedCluster);
            }
        }

        return clusters;
    }

    private static void printHits(List<Hit> hits, String type) {
        System.out.printf("\n%s Hits (%d total):\n", type, hits.size());
        for (Hit hit : hits) {
            System.out.printf("  Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, TDC: %d, ToT: %d, Time: %.2f ns, %sPhi: %.2f rad\n",
                    hit.sector, hit.layer, hit.component, hit.order, hit.tdc, hit.tot, hit.time,
                    (hit.layer != 0 ? String.format("Z: %.2f mm, ", hit.z) : ""), hit.phi);
        }
    }

    private static void printClusters(List<Cluster> clusters, String type) {
        System.out.printf("\n%s Clusters (%d total):\n", type, clusters.size());
        for (Cluster cluster : clusters) {
            System.out.printf("  Cluster ID: %d -> %sPhi: %.2f rad, Time: %.2f ns, Energy: %.2f MeV, Size: %d\n",
                    cluster.id, (cluster.includeZ ? String.format("Z: %.2f mm, ", cluster.z) : ""), cluster.phi, cluster.time, cluster.energy, cluster.hits.size());
            for (Hit hit : cluster.hits) {
                System.out.printf("    Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, TDC: %d, ToT: %d, Time: %.2f ns, %sPhi: %.2f rad\n",
                        hit.sector, hit.layer, hit.component, hit.order, hit.tdc, hit.tot, hit.time,
                        (hit.layer != 0 ? String.format("Z: %.2f mm, ", hit.z) : ""), hit.phi);
            }
        }
    }

    private static void storeClusters(Bank recBank, List<Cluster> clusters) {
        for (int i = 0; i < clusters.size(); i++) {
            Cluster cluster = clusters.get(i);
            recBank.putShort("id", i, (short) cluster.id);
            recBank.putShort("nhits", i, (short) cluster.hits.size());
            recBank.putFloat("phi", i, (float) cluster.phi);
            recBank.putFloat("time", i, (float) cluster.time);
            recBank.putFloat("energy", i, (float) cluster.energy);
            if (cluster.includeZ) {
                recBank.putFloat("z", i, (float) cluster.z);
            }
        }
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
        int id;
        double z, phi, time, energy;
        boolean includeZ;
        List<Hit> hits = new ArrayList<>();

        Cluster(double z, double phi, double time, double energy, boolean includeZ) {
            this.z = z;
            this.phi = phi;
            this.time = time;
            this.energy = energy;
            this.includeZ = includeZ;
        }
    }
}


