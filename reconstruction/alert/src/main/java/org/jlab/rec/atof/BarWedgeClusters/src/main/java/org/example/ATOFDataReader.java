

package org.jlab.rec.atof.BarWedgeClusters;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class ATOFDataReader {

    private static final float VEFF = 20.0f; // Effective speed of light in cm/ns
    private static final float WEDGE_SPACING = 3.0f; // mm
    private static final int WEDGES_PER_BAR = 10;
    private static final int NUM_BARS = 60;
    private static final float L_BAR = 280.0f; // Bar length in mm
    private static final float Z_THRESHOLD = 40.0f; // mm
    private static final float T_THRESHOLD = 5.0f; // ns
    private static final float PHI_THRESHOLD = 0.1f; // radians

    public static void main(String[] args) {
        if (args.length < 3) {
            System.err.println("Usage: java ATOFDataReader <input.hipo> <output.hipo> <schema.json>");
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
        Bank recBank = new Bank(recSchema);

        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(adcBank);

            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();
            extractHits(adcBank, barHits, wedgeHits);

            List<Cluster> barClusters = clusterBarHits(barHits);
            barClusters.forEach(cluster -> System.out.printf("Bar Cluster: Z=%.2f, Phi=%.2f, Hits=%d\n",
                    cluster.z, cluster.phi, cluster.hits.size()));

            List<Cluster> barWedgeClusters = clusterBarWedgeHits(barHits, wedgeHits);
            barWedgeClusters.forEach(cluster -> System.out.printf("Bar+Wedge Cluster: Z=%.2f, Phi=%.2f, Hits=%d\n",
                    cluster.z, cluster.phi, cluster.hits.size()));

            storeClusters(recBank, barClusters, barWedgeClusters);
        }

        reader.close();
        System.out.println("Processing complete!");
    }

    private static void extractHits(Bank adcBank, List<Hit> barHits, List<Hit> wedgeHits) {
        for (int i = 0; i < adcBank.getRows(); i++) {
            int layer = adcBank.getByte("layer", i);
            int component = adcBank.getShort("component", i);
            float time = adcBank.getFloat("time", i);
            int adc = adcBank.getInt("ADC", i);
            if (layer == 0) {
                barHits.add(new Hit(layer, component, time, adc));
            } else if (layer >= 10 && layer <= 19) {
                wedgeHits.add(new Hit(layer, component, time, adc));
            }
        }
    }

    private static List<Cluster> clusterBarHits(List<Hit> barHits) {
        return barHits.stream()
                .collect(Collectors.groupingBy(hit -> hit.component))
                .values().stream()
                .map(Cluster::new)
                .collect(Collectors.toList());
    }

    private static List<Cluster> clusterBarWedgeHits(List<Hit> barHits, List<Hit> wedgeHits) {
        List<Cluster> clusters = new ArrayList<>();
        Map<Integer, List<Hit>> groupedBars = barHits.stream()
                .collect(Collectors.groupingBy(hit -> hit.component));

        for (Hit wedge : wedgeHits) {
            List<Hit> matchingBars = groupedBars.getOrDefault(wedge.component, Collections.emptyList());
            if (matchingBars.size() >= 2) {
                List<Hit> combinedHits = new ArrayList<>(matchingBars);
                combinedHits.add(wedge);
                clusters.add(new Cluster(combinedHits));
            }
        }
        return clusters;
    }

    private static void storeClusters(Bank recBank, List<Cluster> barClusters, List<Cluster> barWedgeClusters) {
        int clusterId = 0;
        for (Cluster cluster : barClusters) {
            storeCluster(recBank, clusterId++, cluster, "Bar");
        }
        for (Cluster cluster : barWedgeClusters) {
            storeCluster(recBank, clusterId++, cluster, "Bar+Wedge");
        }
    }

    private static void storeCluster(Bank recBank, int clusterId, Cluster cluster, String type) {
        recBank.putShort("id", clusterId, (short) clusterId);
        recBank.putShort("nhits", clusterId, (short) cluster.hits.size());
        recBank.putFloat("z", clusterId, cluster.z);
        recBank.putFloat("phi", clusterId, cluster.phi);
        recBank.putFloat("time", clusterId, cluster.minTime);
        recBank.putFloat("energy", clusterId, cluster.energy);
        System.out.printf("Stored Cluster [%s]: ID=%d, Hits=%d, Z=%.2f, Phi=%.2f, Time=%.2f, Energy=%.2f%n",
                type, clusterId, cluster.hits.size(), cluster.z, cluster.phi, cluster.minTime, cluster.energy);
    }

    static class Hit {
        int layer, component, adc;
        float time, z, phi;

        Hit(int layer, int component, float time, int adc) {
            this.layer = layer;
            this.component = component;
            this.time = time;
            this.adc = adc;
            this.z = calculateZ(layer, time);
            this.phi = calculatePhi(component);
        }

        private float calculateZ(int layer, float time) {
            return layer == 0 ? time * VEFF / 2.0f : (component - WEDGES_PER_BAR / 2.0f) * WEDGE_SPACING;
        }

        private float calculatePhi(int component) {
            return (2 * (float) Math.PI * component) / NUM_BARS;
        }
    }

    static class Cluster {
        List<Hit> hits;
        float z, phi, minTime, energy;

        Cluster(List<Hit> hits) {
            this.hits = hits;
            calculateProperties();
        }

        private void calculateProperties() {
            float zSum = 0, phiSum = 0, energySum = 0;
            minTime = Float.MAX_VALUE;

            for (Hit hit : hits) {
                zSum += hit.z;
                phiSum += hit.phi;
                energySum += hit.adc * 0.1f; // Example energy calculation
                minTime = Math.min(minTime, hit.time);
            }

            z = zSum / hits.size();
            phi = phiSum / hits.size();
            energy = energySum;
        }
    }
}









/*
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import java.util.stream.Collectors;

package org.jlab.rec.atof.BarWedgeClusters;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import java.util.*;

public class ATOFDataReader {

    private static final float VEFF = 20.0f; // Effective speed of light in cm/ns
    private static final float WEDGE_SPACING = 3.0f; // mm
    private static final int WEDGES_PER_BAR = 10;
    private static final int NUM_BARS = 60;
    private static final float L_BAR = 280.0f; // Bar length in mm

    private static final float Z_THRESHOLD = 40.0f; // mm
    private static final float T_THRESHOLD = 5.0f; // ns
    private static final float PHI_THRESHOLD = 0.1f; // radians

    public static void main(String[] args) {
        if (args.length < 3) {
            System.err.println("Usage: java ATOFDataReader <input.hipo> <output.hipo> <schema.json>");
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
        Bank recBank = new Bank(recSchema);

        Event event = new Event();

        int eventId = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(adcBank);

            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();
            extractHits(adcBank, barHits, wedgeHits);

            List<Cluster> barClusters = clusterBarHits(barHits);
            for (Cluster cluster : barClusters) {
                System.out.printf("Bar Cluster: Z=%.2f, Phi=%.2f, Hits=%d\n", cluster.z, cluster.phi, cluster.hits.size());
            }

            List<Cluster> barWedgeClusters = clusterBarWedgeHits(barHits, wedgeHits);
            for (Cluster cluster : barWedgeClusters) {
                System.out.printf("Bar+Wedge Cluster: Z=%.2f, Phi=%.2f, Hits=%d\n", cluster.z, cluster.phi, cluster.hits.size());
            }

            storeClusters(recBank, barClusters, barWedgeClusters);

            eventId++;
        }

        reader.close();
        System.out.println("Processing complete!");
    }

    private static void extractHits(Bank adcBank, List<Hit> barHits, List<Hit> wedgeHits) {
        for (int i = 0; i < adcBank.getRows(); i++) {
            int layer = adcBank.getByte("layer", i);
            int component = adcBank.getShort("component", i);
            float time = adcBank.getFloat("time", i);
            int adc = adcBank.getInt("ADC", i);
            if (layer == 0) {
                barHits.add(new Hit(layer, component, time, adc));
            } else if (layer >= 10 && layer <= 19) {
                wedgeHits.add(new Hit(layer, component, time, adc));
            }
        }
    }

    private static List<Cluster> clusterBarHits(List<Hit> barHits) {
        return barHits.stream()
                .collect(Collectors.groupingBy(hit -> hit.component))
                .values().stream()
                .map(Cluster::new)
                .collect(Collectors.toList());
    }

    private static List<Cluster> clusterBarWedgeHits(List<Hit> barHits, List<Hit> wedgeHits) {
        List<Cluster> clusters = new ArrayList<>();
        Map<Integer, List<Hit>> groupedBars = barHits.stream()
                .collect(Collectors.groupingBy(hit -> hit.component));

        for (Hit wedge : wedgeHits) {
            List<Hit> matchingBars = groupedBars.getOrDefault(wedge.component, Collections.emptyList());
            if (matchingBars.size() >= 2) {
                List<Hit> combinedHits = new ArrayList<>(matchingBars);
                combinedHits.add(wedge);
                clusters.add(new Cluster(combinedHits));
            }
        }
        return clusters;
    }

    private static void storeClusters(Bank recBank, List<Cluster> barClusters, List<Cluster> barWedgeClusters) {
        int clusterId = 0;
        for (Cluster cluster : barClusters) {
            storeCluster(recBank, clusterId++, cluster, "Bar");
        }
        for (Cluster cluster : barWedgeClusters) {
            storeCluster(recBank, clusterId++, cluster, "Bar+Wedge");
        }
    }

    private static void storeCluster(Bank recBank, int clusterId, Cluster cluster, String type) {
        recBank.putShort("id", clusterId, (short) clusterId);
        recBank.putShort("nhits", clusterId, (short) cluster.hits.size());
        recBank.putFloat("z", clusterId, cluster.z);
        recBank.putFloat("phi", clusterId, cluster.phi);
        recBank.putFloat("time", clusterId, cluster.minTime);
        recBank.putFloat("energy", clusterId, cluster.energy);
        System.out.printf("Stored Cluster [%s]: ID=%d, Hits=%d, Z=%.2f, Phi=%.2f, Time=%.2f, Energy=%.2f%n",
                type, clusterId, cluster.hits.size(), cluster.z, cluster.phi, cluster.minTime, cluster.energy);
    }

    static class Hit {
        int layer, component, adc;
        float time, z, phi;

        Hit(int layer, int component, float time, int adc) {
            this.layer = layer;
            this.component = component;
            this.time = time;
            this.adc = adc;
            this.z = calculateZ(layer, time);
            this.phi = calculatePhi(component);
        }

        private float calculateZ(int layer, float time) {
            return layer == 0 ? time * VEFF / 2.0f : (component - WEDGES_PER_BAR / 2.0f) * WEDGE_SPACING;
        }

        private float calculatePhi(int component) {
            return (2 * (float) Math.PI * component) / NUM_BARS;
        }
    }

    static class Cluster {
        List<Hit> hits;
        float z, phi, minTime, energy;

        Cluster(List<Hit> hits) {
            this.hits = hits;
            calculateProperties();
        }

        private void calculateProperties() {
            float zSum = 0, phiSum = 0, energySum = 0;
            minTime = Float.MAX_VALUE;

            for (Hit hit : hits) {
                zSum += hit.z;
                phiSum += hit.phi;
                energySum += hit.adc * 0.1f; // Example energy calculation
                minTime = Math.min(minTime, hit.time);
            }

            z = zSum / hits.size();
            phi = phiSum / hits.size();
            energy = energySum;
        }
    }
}
*/
