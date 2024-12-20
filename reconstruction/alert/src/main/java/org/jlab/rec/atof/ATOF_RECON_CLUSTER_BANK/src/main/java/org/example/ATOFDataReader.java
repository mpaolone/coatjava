
package org.example;

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






//works for  with adc
/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.io.HipoWriter;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;

public class ATOFDataReader {

    private static final float DELTA_Z_THRESHOLD = 280.0f; // mm
    private static final float DELTA_PHI_THRESHOLD = 0.1f; // radians
    private static final float DELTA_TIME_THRESHOLD = 2.0f; // ns
    private static final float VEFF = 20.0f; // Effective speed of light in the bar (cm/ns)

    public static void main(String[] args) {
        if (args.length < 3) {
            System.err.println("Usage: java ATOFDataReader <input.hipo> <output.hipo> <json_schema_file>");
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
            System.err.println("Error loading schema from JSON file: " + e.getMessage());
            reader.close();
            System.exit(1);
        }

        if (!schemaFactory.hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the input file.");
            reader.close();
            System.exit(1);
        }

        if (!schemaFactory.hasSchema("ATOF::rec")) {
            System.err.println("Schema ATOF::rec not found. Ensure the schema is loaded from the JSON file.");
            reader.close();
            System.exit(1);
        }

        Schema recSchema = schemaFactory.getSchema("ATOF::rec");
        Bank adcBank = new Bank(schemaFactory.getSchema("ATOF::adc"));
        Event event = new Event();

        HipoWriter writer = new HipoWriter();
        writer.getSchemaFactory().addSchema(recSchema);
        writer.open(outputHipoFile);

        int totalClusters = 0;
        int totalEvents = 0;

        try (PrintWriter logWriter = new PrintWriter(new File("clusters_output.txt"))) {
            while (reader.hasNext()) {
                reader.nextEvent(event);
                event.read(adcBank);

                Map<Integer, List<Hit>> clusters = clusterHits(adcBank);

                int clusterId = 0;
                Bank recBank = new Bank(recSchema, clusters.size());

                for (Map.Entry<Integer, List<Hit>> entry : clusters.entrySet()) {
                    List<Hit> clusterHits = entry.getValue();

                    boolean hasBarHits = clusterHits.stream().anyMatch(hit -> hit.layer == 0);
                    boolean hasWedgeHits = clusterHits.stream().anyMatch(hit -> hit.layer >= 10 && hit.layer <= 19);

                    if (hasBarHits && hasWedgeHits && clusterHits.size() < 3) continue;
                    if (hasBarHits && !hasWedgeHits && clusterHits.size() < 2) continue;

                    int totalADC = 0;
                    float weightedPhi = 0;
                    float minTime = Float.MAX_VALUE;
                    float totalEnergy = 0;
                    float zCluster = 0;

                    float tLeft = Float.MAX_VALUE;
                    float tRight = Float.MAX_VALUE;

                    for (Hit hit : clusterHits) {
                        totalADC += hit.adc;
                        weightedPhi += hit.phi * hit.adc;
                        totalEnergy += hit.adc * 0.1f;
                        minTime = Math.min(minTime, hit.time);

                        if (hit.layer == 0) {
                            if (hit.component % 2 == 0) {
                                tLeft = Math.min(tLeft, hit.time);
                            } else {
                                tRight = Math.min(tRight, hit.time);
                            }
                        }
                    }

                    if (totalADC > 0) {
                        weightedPhi /= totalADC;
                    } else {
                        weightedPhi = Float.NaN;
                    }

                    if (hasBarHits) {
                        // ZBar dominates the cluster Z position
                        if (tLeft != Float.MAX_VALUE && tRight != Float.MAX_VALUE) {
                            zCluster = VEFF * (tLeft - tRight) / 2.0f;
                        } else if (tLeft != Float.MAX_VALUE) {
                            zCluster = VEFF * tLeft;
                        } else if (tRight != Float.MAX_VALUE) {
                            zCluster = VEFF * tRight;
                        } else {
                            zCluster = Float.NaN;
                        }
                    } else {
                        // For wedge-only clusters, average Z of wedges
                        zCluster = (float) clusterHits.stream()
                                .filter(hit -> hit.layer >= 10 && hit.layer <= 19)
                                .mapToDouble(hit -> hit.z)
                                .average()
                                .orElse(Float.NaN);
                    }

                    recBank.putShort("id", clusterId, (short) clusterId);
                    recBank.putShort("nhits", clusterId, (short) clusterHits.size());
                    recBank.putFloat("z", clusterId, zCluster);
                    recBank.putFloat("phi", clusterId, weightedPhi);
                    recBank.putFloat("time", clusterId, minTime);
                    recBank.putFloat("energy", clusterId, totalEnergy);

                    clusterId++;
                }

                totalClusters += recBank.getRows();
                totalEvents++;

                logWriter.printf("Event %d: %d Clusters%n", totalEvents, recBank.getRows());
                for (int row = 0; row < recBank.getRows(); row++) {
                    logWriter.printf("Cluster ID: %d, Hits: %d, Z: %.3f mm, Phi: %.3f rad, Time: %.3f ns, Energy: %.3f MeV%n",
                            recBank.getShort("id", row),
                            recBank.getShort("nhits", row),
                            recBank.getFloat("z", row),
                            recBank.getFloat("phi", row),
                            recBank.getFloat("time", row),
                            recBank.getFloat("energy", row));
                }

                event.write(recBank);
                writer.addEvent(event);
            }

            logWriter.printf("Total Events: %d%n", totalEvents);
            logWriter.printf("Total Clusters: %d%n", totalClusters);
        } catch (Exception e) {
            System.err.println("Error writing log file: " + e.getMessage());
        }

        reader.close();
        writer.close();

        System.out.printf("Processing complete. Total events: %d, Total clusters: %d. Output written to: %s%n", totalEvents, totalClusters, outputHipoFile);
    }

    private static Map<Integer, List<Hit>> clusterHits(Bank adcBank) {
        Map<Integer, List<Hit>> clusters = new HashMap<>();
        List<Hit> hits = new ArrayList<>();

        for (int i = 0; i < adcBank.getRows(); i++) {
            int layer = adcBank.getByte("layer", i);
            int adc = adcBank.getInt("ADC", i);
            float time = adcBank.getFloat("time", i);
            int component = adcBank.getShort("component", i);

            float z = calculateZ(layer, component);
            float phi = calculatePhi(component);

            hits.add(new Hit(i, z, phi, time, adc, layer, component));
        }

        for (Hit hit : hits) {
            boolean addedToCluster = false;
            for (List<Hit> cluster : clusters.values()) {
                if (cluster.stream().anyMatch(existing -> isWithinThreshold(existing, hit))) {
                    cluster.add(hit);
                    addedToCluster = true;
                    break;
                }
            }
            if (!addedToCluster) {
                clusters.put(clusters.size(), new ArrayList<>(Collections.singletonList(hit)));
            }
        }

        return clusters;
    }

    private static boolean isWithinThreshold(Hit h1, Hit h2) {
        return Math.abs(h1.z - h2.z) < DELTA_Z_THRESHOLD &&
                Math.abs(h1.phi - h2.phi) < DELTA_PHI_THRESHOLD &&
                Math.abs(h1.time - h2.time) < DELTA_TIME_THRESHOLD;
    }

    private static float calculateZ(int layer, int component) {
        if (layer == 0) return 0.0f;
        return (component % 10 - 4.5f) * 3.0f;
    }

    private static float calculatePhi(int component) {
        return (float) (2.0 * Math.PI * component / 60.0);
    }

    static class Hit {
        int id, adc, layer, component;
        float z, phi, time;

        Hit(int id, float z, float phi, float time, int adc, int layer, int component) {
            this.id = id;
            this.z = z;
            this.phi = phi;
            this.time = time;
            this.adc = adc;
            this.layer = layer;
            this.component = component;
        }
    }
}


*/

