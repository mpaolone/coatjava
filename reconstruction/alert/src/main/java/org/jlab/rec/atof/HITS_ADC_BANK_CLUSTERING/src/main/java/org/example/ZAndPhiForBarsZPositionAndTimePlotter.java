package org.jlab.rec.atof.HITS_ADC_BANK_CLUSTERING;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotter {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm
    private static final double PHI_THRESHOLD = 0.1; // rad
    private static final double TIME_THRESHOLD = 1.0; // ns

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processHitsGlobally(reader);
        reader.close();

        // Print summary of cluster sizes
        System.out.println("\nCluster Size Summary:");
        for (Map.Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            System.out.printf("Clusters of size %d: %d\n", entry.getKey(), entry.getValue());
        }
    }

    private static void processHitsGlobally(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();

        List<Hit> barHits = new ArrayList<>();
        List<Hit> wedgeHits = new ArrayList<>();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofAdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit); // Bar hits
                else if (hit.layer >= 10 && hit.layer <= 19) wedgeHits.add(hit); // Wedge hits
            }
        }

        for (int i = 0; i < barHits.size(); i++) {
            Hit barHit1 = barHits.get(i);

            for (int j = i + 1; j < barHits.size(); j++) {
                Hit barHit2 = barHits.get(j);

                if (barHit1.sector == barHit2.sector &&
                    barHit1.layer == barHit2.layer &&
                    barHit1.component == barHit2.component) {

                    double deltaT = barHit1.time - barHit2.time;
                    double zBar = VEFF * deltaT / 2;
                    double tBar = Math.min(barHit1.time - (zBar - BAR_LENGTH / 2) / VEFF,
                                           barHit2.time - (zBar + BAR_LENGTH / 2) / VEFF);

                    List<Hit> clusterWedgeHits = new ArrayList<>();
                    int clusterSize = 2;

                    for (Hit wedgeHit : wedgeHits) {
                        double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                        double deltaPhi = Math.abs(barHit1.phi - wedgeHit.phi);
                        double deltaTime = Math.abs(tBar - wedgeHit.time);

                        deltaZList.add(deltaZ);
                        deltaPhiList.add(deltaPhi);
                        deltaTimeList.add(deltaTime);

                        if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                            clusterSize++;
                            clusterWedgeHits.add(wedgeHit);
                        }
                    }

                    clusterSizes.add(clusterSize);
                    clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

                    // Print cluster information
                    System.out.printf("Cluster Formed (Size %d):\n", clusterSize);
                    System.out.printf("  Bar Hits:\n");
                    printHitDetails(barHit1, "Bar Hit 1");
                    printHitDetails(barHit2, "Bar Hit 2");
                    System.out.printf("  Combined Bar Z and T: ZBar: %.2f mm, TBar: %.2f ns\n", zBar, tBar);

                    System.out.println("  Wedge Hits:");
                    for (Hit wedgeHit : clusterWedgeHits) {
                        printHitDetails(wedgeHit, "Wedge Hit");
                    }
                    System.out.println("---- End of Cluster ----\n");
                }
            }
        }
    }

    private static Hit createHit(Bank bank, int index) {
        int sector = bank.getByte("sector", index);
        int layer = bank.getByte("layer", index);
        int component = bank.getShort("component", index);
        int order = bank.getByte("order", index);
        int adc = bank.getInt("ADC", index);
        float time = bank.getFloat("time", index);
        double phi = 2 * Math.PI * component / 60;
        return new Hit(sector, layer, component, order, adc, time, phi);
    }

    private static void printHitDetails(Hit hit, String label) {
        System.out.printf("%s -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad, Z: %.2f mm)\n",
                label, hit.sector, hit.layer, hit.component, hit.order, hit.adc, hit.time, hit.phi, hit.zWedge());
    }

    static class Hit {
        int sector, layer, component, order, adc;
        double time, phi;

        Hit(int sector, int layer, int component, int order, int adc, double time, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.time = time;
            this.phi = phi;
        }

        double zWedge() {
            int wedgeIndex = component % N_WEDGE;
            return (wedgeIndex - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;
        }
    }
}


/*

package org.jlab.rec.atof.HITS_TDC_BANK_CLUSTERING;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotter {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm
    private static final double PHI_THRESHOLD = 0.1; // rad
    private static final double TIME_THRESHOLD = 1.0; // ns
    private static final double TDC_RESOLUTION = 0.015625; // ns per TDC unit

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::tdc")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processHitsGlobally(reader);
        reader.close();

        // Print summary of cluster sizes
        System.out.println("\nCluster Size Summary:");
        for (Map.Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            System.out.printf("Clusters of size %d: %d\n", entry.getKey(), entry.getValue());
        }
    }

    private static void processHitsGlobally(HipoReader reader) {
        Bank atofTdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::tdc"));
        Event event = new Event();

        List<Hit> barHits = new ArrayList<>();
        List<Hit> wedgeHits = new ArrayList<>();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofTdcBank);

            int numHits = atofTdcBank.getRows();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofTdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit); // Bar hits
                else if (hit.layer >0 && hit.layer <4) wedgeHits.add(hit); // Wedge hits
            }
        }

        for (int i = 0; i < barHits.size(); i++) {
            Hit barHit1 = barHits.get(i);

            for (int j = i + 1; j < barHits.size(); j++) {
                Hit barHit2 = barHits.get(j);

                if (barHit1.sector == barHit2.sector &&
                    barHit1.layer == barHit2.layer &&
                    barHit1.component == barHit2.component) {

                    double deltaT = barHit1.time - barHit2.time;
                    double zBar = VEFF * deltaT / 2;
                    double tBar = Math.min(barHit1.time - (zBar - BAR_LENGTH / 2) / VEFF,
                                           barHit2.time - (zBar + BAR_LENGTH / 2) / VEFF);

                    List<Hit> clusterWedgeHits = new ArrayList<>();
                    int clusterSize = 2;

                    for (Hit wedgeHit : wedgeHits) {
                        double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                        double deltaPhi = Math.abs(barHit1.phi - wedgeHit.phi);
                        double deltaTime = Math.abs(tBar - wedgeHit.time);

                        deltaZList.add(deltaZ);
                        deltaPhiList.add(deltaPhi);
                        deltaTimeList.add(deltaTime);

                        if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                            clusterSize++;
                            clusterWedgeHits.add(wedgeHit);
                        }
                    }

                    clusterSizes.add(clusterSize);
                    clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

                    // Print cluster information
                    System.out.printf("Cluster Formed (Size %d):\n", clusterSize);
                    System.out.printf("  Bar Hits:\n");
                    printHitDetails(barHit1, "Bar Hit 1");
                    printHitDetails(barHit2, "Bar Hit 2");
                    System.out.printf("  Combined Bar Z and T: ZBar: %.2f mm, TBar: %.2f ns\n", zBar, tBar);

                    System.out.println("  Wedge Hits:");
                    for (Hit wedgeHit : clusterWedgeHits) {
                        printHitDetails(wedgeHit, "Wedge Hit");
                    }
                    System.out.println("---- End of Cluster ----\n");
                }
            }
        }
    }

    private static Hit createHit(Bank bank, int index) {
        int sector = bank.getByte("sector", index);
        int layer = bank.getByte("layer", index);
        int component = bank.getShort("component", index);
        int order = bank.getByte("order", index);
        int tdc = bank.getInt("TDC", index);
        double time = tdc * TDC_RESOLUTION;
        double phi = 2 * Math.PI * component / 60;
        return new Hit(sector, layer, component, order, time, phi);
    }

    private static void printHitDetails(Hit hit, String label) {
        System.out.printf("%s -> (Sector: %d, Layer: %d, Component: %d, Order: %d, Time: %.2f ns, Phi: %.2f rad, Z: %.2f mm)\n",
                label, hit.sector, hit.layer, hit.component, hit.order, hit.time, hit.phi, hit.zWedge());
    }

    static class Hit {
        int sector, layer, component, order;
        double time, phi;

        Hit(int sector, int layer, int component, int order, double time, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.time = time;
            this.phi = phi;
        }

        double zWedge() {
            int wedgeIndex = component % N_WEDGE;
            return (wedgeIndex - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;
        }
    }
}
*/
