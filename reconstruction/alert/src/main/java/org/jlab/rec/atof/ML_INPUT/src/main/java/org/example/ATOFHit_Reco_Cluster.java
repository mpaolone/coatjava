package org.jlab.rec.atof.ML_INPUT;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class ATOFHit_Reco_Cluster {


    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10;
    private static final double Z_THRESHOLD = 280.0;
    private static final double PHI_THRESHOLD = 0.1; // radians
    private static final double TIME_THRESHOLD = 1.7;


    private static List<String[]> csvData = new ArrayList<>();

    public static void main(String[] args) {
        if (args.length < 2) {
            System.err.println("Usage: java ATOFHit_Reco_Cluster <hipoFilePath> <outputCSVFile>");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        String csvFilePath = args[1];

        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc") || !reader.getSchemaFactory().hasSchema("MC::True")) {
            System.err.println("Required schemas not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        csvData.add(new String[]{"EventID", "HitID", "ClusterID", "DeltaZ", "DeltaPhi", "DeltaTime", "ADC", "ZPosition", "Phi", "Time"});

        processEvents(reader);
        reader.close();


        writeToCSV(csvFilePath);

        System.out.println("CSV file written to: " + csvFilePath);
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Bank mcTrueBank = new Bank(reader.getSchemaFactory().getSchema("MC::True"));
        Event event = new Event();

        int eventCount = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);
            event.read(mcTrueBank);

            int numHits = atofAdcBank.getRows();
            List<Hit> hits = new ArrayList<>();

            System.out.printf("\nProcessing Event #%d with %d hits.\n", eventCount, numHits);

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                hits.add(createHit(atofAdcBank, hitIndex));
            }


            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();
            for (Hit hit : hits) {
                if (hit.getLayer() == 0) barHits.add(hit);
                else if (hit.getLayer() >= 10 && hit.getLayer() <= 19) wedgeHits.add(hit);
            }

            if (barHits.isEmpty() && wedgeHits.isEmpty()) {
                System.out.println("No valid bar or wedge hits in this event.");
                eventCount++;
                continue;
            }

            if (barHits.size() >= 2) {
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);

                double zBar = VEFF * (barLeft.getTime() - barRight.getTime()) / 2;
                double tBar = Math.min(barLeft.getTime() - (zBar - BAR_LENGTH / 2) / VEFF,
                        barRight.getTime() - (zBar + BAR_LENGTH / 2) / VEFF);

                processCluster(barLeft, barRight, wedgeHits, zBar, tBar, eventCount);
            }

            eventCount++;
        }
    }

    private static void processCluster(Hit barLeft, Hit barRight, List<Hit> wedgeHits, double zBar, double tBar, int eventID) {
        List<Hit> clusterWedgeHits = new ArrayList<>();
        int clusterID = 1;

        for (Hit wedgeHit : wedgeHits) {
            double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
            double deltaPhi = normalizePhi(barLeft.getPhi() - wedgeHit.getPhi());
            double deltaTime = Math.abs(tBar - wedgeHit.getTime());

            if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                clusterWedgeHits.add(wedgeHit);
                csvData.add(new String[]{
                        String.valueOf(eventID),
                        String.valueOf(wedgeHit.getComponent()),
                        String.valueOf(clusterID),
                        String.format("%.2f", deltaZ),
                        String.format("%.2f", deltaPhi),
                        String.format("%.2f", deltaTime),
                        String.valueOf(wedgeHit.getAdc()),
                        String.format("%.2f", wedgeHit.zWedge()),
                        String.format("%.2f", wedgeHit.getPhi()),
                        String.format("%.2f", wedgeHit.getTime())
                });
            }
        }


        csvData.add(new String[]{
                String.valueOf(eventID),
                String.valueOf(barLeft.getComponent()),
                String.valueOf(clusterID),
                "0.0",
                "0.0",
                "0.0",
                String.valueOf(barLeft.getAdc()),
                String.format("%.2f", zBar),
                String.format("%.2f", barLeft.getPhi()),
                String.format("%.2f", barLeft.getTime())
        });

        csvData.add(new String[]{
                String.valueOf(eventID),
                String.valueOf(barRight.getComponent()),
                String.valueOf(clusterID),
                "0.0",
                "0.0",
                "0.0",
                String.valueOf(barRight.getAdc()),
                String.format("%.2f", zBar),
                String.format("%.2f", barRight.getPhi()),
                String.format("%.2f", barRight.getTime())
        });
    }

    private static Hit createHit(Bank bank, int index) {
        int sector = bank.getByte("sector", index);
        int layer = bank.getByte("layer", index);
        int component = bank.getShort("component", index);
        int order = bank.getByte("order", index);
        int adc = bank.getInt("ADC", index);
        int ped = bank.getInt("ped", index);
        float time = bank.getFloat("time", index);
        double phi = 2 * Math.PI * component / 60;
        return new Hit(sector, layer, component, order, adc, time, phi, ped);
    }

    private static double normalizePhi(double phi) {
        while (phi > Math.PI) phi -= 2 * Math.PI;
        while (phi < -Math.PI) phi += 2 * Math.PI;
        return phi;
    }

    private static void writeToCSV(String filePath) {
        try (FileWriter writer = new FileWriter(filePath)) {
            for (String[] line : csvData) {
                writer.write(String.join(",", line));
                writer.write("\n");
            }
        } catch (IOException e) {
            System.err.println("Error writing to CSV file: " + e.getMessage());
        }
    }

    static class Hit {
        private int sector, layer, component, order, adc, ped;
        private double time, phi;

        Hit(int sector, int layer, int component, int order, int adc, double time, double phi, int ped) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.time = time;
            this.phi = phi;
            this.ped = ped;
        }

        public int getSector() { return sector; }
        public int getLayer() { return layer; }
        public int getComponent() { return component; }
        public int getOrder() { return order; }
        public int getAdc() { return adc; }
        public double getTime() { return time; }
        public double getPhi() { return phi; }
        public int getPed() { return ped; }

        double zWedge() {
            int wedgeIndex = component % N_WEDGE;
            return (wedgeIndex - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;
        }
    }
}

