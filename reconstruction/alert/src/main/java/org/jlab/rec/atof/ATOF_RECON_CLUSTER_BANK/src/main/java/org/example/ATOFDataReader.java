package org.jlab.rec.atof.ATOF_RECON_CLUSTER_BANK;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.io.HipoWriter;

import java.util.*;

public class ATOFDataReader {
    private static final float DELTA_Z_THRESHOLD = 280.0f; // mm
    private static final float DELTA_PHI_THRESHOLD = 0.1f; // radians
    private static final float DELTA_TIME_THRESHOLD = 2.0f; // ns
    private static final float VEFF = 20.0f; // Effective speed of light in the bar (cm/ns)

    public static void main(String[] args) {
        if (args.length < 3) {
            System.err.println("Usage: java ATOFClusterReconstruction <input.hipo> <output.hipo> <json_schema_file>");
            System.exit(1);
        }

        String inputHipoFile = args[0];
        String outputHipoFile = args[1];
        String jsonSchemaFile = args[2];

        HipoReader reader = new HipoReader();
        reader.open(inputHipoFile);

        SchemaFactory schemaFactory = reader.getSchemaFactory();
        try {
            schemaFactory.readFile(jsonSchemaFile);
        } catch (Exception e) {
            System.err.println("Error loading schema: " + e.getMessage());
            reader.close();
            System.exit(1);
        }

        Schema recSchema = schemaFactory.getSchema("ATOF::rec");
        Bank adcBank = new Bank(schemaFactory.getSchema("ATOF::adc"));
        Event event = new Event();

        HipoWriter writer = new HipoWriter();
        writer.getSchemaFactory().addSchema(recSchema);
        writer.open(outputHipoFile);

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(adcBank);

            List<Cluster> clusters = reconstructClusters(adcBank);

            Bank recBank = new Bank(recSchema, clusters.size());
            for (int i = 0; i < clusters.size(); i++) {
                Cluster cluster = clusters.get(i);

                recBank.putShort("id", i, (short) i);
                recBank.putShort("nhits", i, (short) cluster.hits.size());
                recBank.putFloat("z", i, cluster.z);
                recBank.putFloat("phi", i, cluster.phi);
                recBank.putFloat("time", i, cluster.time);
                recBank.putFloat("energy", i, cluster.energy);
                recBank.putByte("sector", i, (byte) cluster.sector);
                recBank.putByte("layer", i, (byte) cluster.layer);
                recBank.putShort("component", i, (short) cluster.component);
                recBank.putByte("order", i, (byte) cluster.order);
                recBank.putFloat("ADC", i, cluster.adc);
                recBank.putFloat("ped", i, cluster.ped);
                recBank.putFloat("ZBar", i, cluster.zBar);
                recBank.putFloat("ZWedge", i, cluster.zWedge);
                recBank.putFloat("deltaZ", i, cluster.deltaZ);
                recBank.putFloat("deltaPhi", i, cluster.deltaPhi);
                recBank.putFloat("deltaTime", i, cluster.deltaTime);
                recBank.putShort("barHits", i, (short) cluster.barHits);
                recBank.putShort("wedgeHits", i, (short) cluster.wedgeHits);
            }

            event.write(recBank);
            writer.addEvent(event);
        }

        reader.close();
        writer.close();

        System.out.printf("Processing complete. Output written to: %s%n", outputHipoFile);
    }

    private static List<Cluster> reconstructClusters(Bank adcBank) {
        List<Cluster> clusters = new ArrayList<>();
        List<Hit> hits = new ArrayList<>();

        for (int i = 0; i < adcBank.getRows(); i++) {
            int sector = adcBank.getByte("sector", i);
            int layer = adcBank.getByte("layer", i);
            int component = adcBank.getShort("component", i);
            int order = adcBank.getByte("order", i);
            float adc = adcBank.getFloat("ADC", i);
            float ped = adcBank.getFloat("ped", i);
            float time = adcBank.getFloat("time", i);

            float z = calculateZ(layer, component, time, order);
            float phi = calculatePhi(component);

            hits.add(new Hit(sector, layer, component, order, adc, ped, time, z, phi));
        }

        for (Hit hit : hits) {
            boolean addedToCluster = false;
            for (Cluster cluster : clusters) {
                if (cluster.canAddHit(hit)) {
                    cluster.addHit(hit);
                    addedToCluster = true;
                    break;
                }
            }

            if (!addedToCluster) {
                clusters.add(new Cluster(hit));
            }
        }

        return clusters;
    }

    private static float calculateZ(int layer, int component, float time, int order) {
        if (layer == 0) return VEFF * (order == 0 ? time : -time);
        return (component % 10 - 4.5f) * 3.0f; // Wedge Z
    }

    private static float calculatePhi(int component) {
        return (float) (2.0 * Math.PI * component / 60.0);
    }

    static class Hit {
        int sector, layer, component, order;
        float adc, ped, time, z, phi;

        Hit(int sector, int layer, int component, int order, float adc, float ped, float time, float z, float phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.ped = ped;
            this.time = time;
            this.z = z;
            this.phi = phi;
        }
    }

    static class Cluster {
        List<Hit> hits = new ArrayList<>();
        float z, phi, time, adc, ped, energy, zBar, zWedge, deltaZ, deltaPhi, deltaTime;
        int sector, layer, component, order;
        int barHits = 0, wedgeHits = 0;

        Cluster(Hit hit) {
            addHit(hit);
        }

        void addHit(Hit hit) {
            hits.add(hit);
            adc += hit.adc;
            ped += hit.ped;
            energy += hit.adc * 0.1f;
            z += hit.z;
            phi += hit.phi;
            time = Math.min(time, hit.time);
            if (hit.layer == 0) barHits++;
            else wedgeHits++;
        }

        boolean canAddHit(Hit hit) {
            return Math.abs(hit.z - z) < DELTA_Z_THRESHOLD &&
                    Math.abs(hit.phi - phi) < DELTA_PHI_THRESHOLD &&
                    Math.abs(hit.time - time) < DELTA_TIME_THRESHOLD;
        }
    }
}


/*

//equivalent TDC version of above code: 

package org.jlab.rec.atof.ATOF_RECON_CLUSTER_BANK;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.io.HipoWriter;

import java.util.*;

public class ATOFDataReader {
    private static final float DELTA_Z_THRESHOLD = 280.0f; // mm
    private static final float DELTA_PHI_THRESHOLD = 0.1f; // radians
    private static final float DELTA_TIME_THRESHOLD = 2.0f; // ns
    private static final float VEFF = 20.0f; // Effective speed of light in the bar (cm/ns)
    private static final float TDC_RESOLUTION = 0.015625f; // TDC resolution in ns

    public static void main(String[] args) {
        if (args.length < 3) {
            System.err.println("Usage: java ATOFClusterReconstruction <input.hipo> <output.hipo> <json_schema_file>");
            System.exit(1);
        }

        String inputHipoFile = args[0];
        String outputHipoFile = args[1];
        String jsonSchemaFile = args[2];

        HipoReader reader = new HipoReader();
        reader.open(inputHipoFile);

        SchemaFactory schemaFactory = reader.getSchemaFactory();
        try {
            schemaFactory.readFile(jsonSchemaFile);
        } catch (Exception e) {
            System.err.println("Error loading schema: " + e.getMessage());
            reader.close();
            System.exit(1);
        }

        Schema recSchema = schemaFactory.getSchema("ATOF::rec");
        Bank tdcBank = new Bank(schemaFactory.getSchema("ATOF::tdc"));
        Event event = new Event();

        HipoWriter writer = new HipoWriter();
        writer.getSchemaFactory().addSchema(recSchema);
        writer.open(outputHipoFile);

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(tdcBank);

            List<Cluster> clusters = reconstructClusters(tdcBank);

            Bank recBank = new Bank(recSchema, clusters.size());
            for (int i = 0; i < clusters.size(); i++) {
                Cluster cluster = clusters.get(i);

                recBank.putShort("id", i, (short) i);
                recBank.putShort("nhits", i, (short) cluster.hits.size());
                recBank.putFloat("z", i, cluster.z);
                recBank.putFloat("phi", i, cluster.phi);
                recBank.putFloat("time", i, cluster.time);
                recBank.putFloat("energy", i, cluster.energy);
                recBank.putByte("sector", i, (byte) cluster.sector);
                recBank.putByte("layer", i, (byte) cluster.layer);
                recBank.putShort("component", i, (short) cluster.component);
                recBank.putByte("order", i, (byte) cluster.order);
                recBank.putFloat("TDC", i, cluster.tdc);
                recBank.putFloat("ZBar", i, cluster.zBar);
                recBank.putFloat("ZWedge", i, cluster.zWedge);
                recBank.putFloat("deltaZ", i, cluster.deltaZ);
                recBank.putFloat("deltaPhi", i, cluster.deltaPhi);
                recBank.putFloat("deltaTime", i, cluster.deltaTime);
                recBank.putShort("barHits", i, (short) cluster.barHits);
                recBank.putShort("wedgeHits", i, (short) cluster.wedgeHits);
            }

            event.write(recBank);
            writer.addEvent(event);
        }

        reader.close();
        writer.close();

        System.out.printf("Processing complete. Output written to: %s%n", outputHipoFile);
    }

    private static List<Cluster> reconstructClusters(Bank tdcBank) {
        List<Cluster> clusters = new ArrayList<>();
        List<Hit> hits = new ArrayList<>();

        for (int i = 0; i < tdcBank.getRows(); i++) {
            int sector = tdcBank.getInt("sector", i);
            int layer = tdcBank.getInt("layer", i);
            int component = tdcBank.getInt("component", i);
            int order = tdcBank.getInt("order", i);
            float tdc = tdcBank.getInt("TDC", i) * TDC_RESOLUTION;

            float z = calculateZ(layer, component, tdc, order);
            float phi = calculatePhi(component);

            hits.add(new Hit(sector, layer, component, order, tdc, z, phi));
        }

        for (Hit hit : hits) {
            boolean addedToCluster = false;
            for (Cluster cluster : clusters) {
                if (cluster.canAddHit(hit)) {
                    cluster.addHit(hit);
                    addedToCluster = true;
                    break;
                }
            }

            if (!addedToCluster) {
                clusters.add(new Cluster(hit));
            }
        }

        return clusters;
    }

    private static float calculateZ(int layer, int component, float time, int order) {
        if (layer == 0) return VEFF * (order == 0 ? time : -time);
        return (component % 10 - 4.5f) * 3.0f; // Wedge Z
    }

    private static float calculatePhi(int component) {
        return (float) (2.0 * Math.PI * component / 60.0);
    }

    static class Hit {
        int sector, layer, component, order;
        float tdc, z, phi;

        Hit(int sector, int layer, int component, int order, float tdc, float z, float phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.tdc = tdc;
            this.z = z;
            this.phi = phi;
        }
    }

    static class Cluster {
        List<Hit> hits = new ArrayList<>();
        float z, phi, time, tdc, energy, zBar, zWedge, deltaZ, deltaPhi, deltaTime;
        int sector, layer, component, order;
        int barHits = 0, wedgeHits = 0;

        Cluster(Hit hit) {
            addHit(hit);
        }

        void addHit(Hit hit) {
            hits.add(hit);
            tdc += hit.tdc;
            energy += hit.tdc * 0.1f;
            z += hit.z;
            phi += hit.phi;
            time = Math.min(time, hit.tdc);
            if (hit.layer == 0) barHits++;
            else wedgeHits++;
        }

        boolean canAddHit(Hit hit) {
            return Math.abs(hit.z - z) < DELTA_Z_THRESHOLD &&
                    Math.abs(hit.phi - phi) < DELTA_PHI_THRESHOLD &&
                    Math.abs(hit.tdc - time) < DELTA_TIME_THRESHOLD;
        }
    }
}

*/
