package org.jlab.rec.atof.BAR_WEDGE_CLUSTERING;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import java.util.*;

public class ATOFHit_Reco_Cluster {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 3.0; // mm
    private static final int N_WEDGE = 10;
    private static final double Z_THRESHOLD = 280.0;
    private static final double PHI_THRESHOLD = 0.3;
    private static final double TIME_THRESHOLD = 1.7;

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

        processEvents(reader);
        reader.close();

        System.out.println("\nCluster Size Summary:");
        for (Map.Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            System.out.printf("Clusters of size %d: %d\n", entry.getKey(), entry.getValue());
        }
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();

        int eventCount = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

            System.out.printf("Processing Event #%d with %d hits.\n", eventCount, numHits);

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofAdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit);
                else if (hit.layer >= 10 && hit.layer <= 19) wedgeHits.add(hit);
            }

            if (barHits.size() >= 2) {
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);

                double zBar = VEFF * (barLeft.time - barRight.time) / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF, barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);

                System.out.printf("Bar Hits (for Cluster Calculation):\n");
                printHitDetails(barLeft, "Bar Hit #1");
                printHitDetails(barRight, "Bar Hit #2");
                System.out.printf("  Calculated ZBar: %.2f mm, TBar: %.2f ns\n", zBar, tBar);

                List<Hit> clusterWedgeHits = new ArrayList<>();
                for (Hit wedgeHit : wedgeHits) {
                    double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                    double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        clusterWedgeHits.add(wedgeHit);
                    }
                }

                if (!clusterWedgeHits.isEmpty()) {
                    int clusterSize = 2 + clusterWedgeHits.size();
                    clusterSizes.add(clusterSize);
                    clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

                    System.out.printf("Cluster Formed (Size: %d)\n", clusterSize);
                }
            }
            eventCount++;
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
        System.out.printf("%s -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad, Z: %.2f mm\n",
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

//TDC version 

package org.jlab.rec.atof.BAR_WEDGE_CLUSTERING;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import java.util.*;

public class ATOFHitRecoCluster {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 3.0; // mm
    private static final int N_WEDGE = 10;
    private static final double Z_THRESHOLD = 280.0;
    private static final double PHI_THRESHOLD = 0.3;
    private static final double TIME_THRESHOLD = 1.7;

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
            System.err.println("Required schema ATOF::tdc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        System.out.println("\nCluster Size Summary:");
        for (Map.Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            System.out.printf("Clusters of size %d: %d\n", entry.getKey(), entry.getValue());
        }
    }

    private static void processEvents(HipoReader reader) {
        Bank atofTdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::tdc"));
        Event event = new Event();

        int eventCount = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofTdcBank);

            int numHits = atofTdcBank.getRows();
            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

            System.out.printf("Processing Event #%d with %d hits.\n", eventCount, numHits);

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofTdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit);
                else if (hit.layer >= 10 && hit.layer <= 19) wedgeHits.add(hit);
            }

            if (barHits.size() >= 2) {
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);

                double zBar = VEFF * (barLeft.time - barRight.time) / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF, barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);

                System.out.printf("Bar Hits (for Cluster Calculation):\n");
                printHitDetails(barLeft, "Bar Hit #1");
                printHitDetails(barRight, "Bar Hit #2");
                System.out.printf("  Calculated ZBar: %.2f mm, TBar: %.2f ns\n", zBar, tBar);

                List<Hit> clusterWedgeHits = new ArrayList<>();
                for (Hit wedgeHit : wedgeHits) {
                    double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                    double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        clusterWedgeHits.add(wedgeHit);
                    }
                }

                if (!clusterWedgeHits.isEmpty()) {
                    int clusterSize = 2 + clusterWedgeHits.size();
                    clusterSizes.add(clusterSize);
                    clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

                    System.out.printf("Cluster Formed (Size: %d)\n", clusterSize);
                }
            }
            eventCount++;
        }
    }

    private static Hit createHit(Bank bank, int index) {
        int sector = bank.getByte("sector", index);
        int layer = bank.getByte("layer", index);
        int component = bank.getShort("component", index);
        int order = bank.getByte("order", index);
        int tdc = bank.getInt("TDC", index);
        int tot = bank.getInt("ToT", index);
        double time = tdc * 0.015625; // Convert TDC to ns
        double phi = 2 * Math.PI * component / 60;
        return new Hit(sector, layer, component, order, tdc, tot, time, phi);
    }

    private static void printHitDetails(Hit hit, String label) {
        System.out.printf("%s -> Sector: %d, Layer: %d, Component: %d, Order: %d, TDC: %d, ToT: %d, Time: %.2f ns, Phi: %.2f rad, Z: %.2f mm\n",
                label, hit.sector, hit.layer, hit.component, hit.order, hit.tdc, hit.tot, hit.time, hit.phi, hit.zWedge());
    }

    static class Hit {
        int sector, layer, component, order, tdc, tot;
        double time, phi;

        Hit(int sector, int layer, int component, int order, int tdc, int tot, double time, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.tdc = tdc;
            this.tot = tot;
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
