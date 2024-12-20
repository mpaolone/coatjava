package org.jlab.rec.atof.ML_INPUT.org.example;

//package org.example;

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

/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.statistics.HistogramDataset;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ATOFHit_Reco_Cluster extends JFrame {

    // Constants for detector and clustering thresholds
    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10;
    private static final double Z_THRESHOLD = 280.0;
    private static final double PHI_THRESHOLD = 0.1; // radians
    private static final double TIME_THRESHOLD = 1.7;

    // Lists for reconstructed clustering
    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

    // Lists for cluster properties
    private static List<Double> zClusterList = new ArrayList<>();
    private static List<Double> phiClusterList = new ArrayList<>();
    private static List<Double> timeClusterList = new ArrayList<>();
    private static List<Double> energyClusterList = new ArrayList<>();

    // Lists for residuals before clustering
    private static List<Double> deltaZResiduals_BeforeClustering = new ArrayList<>();
    private static List<Double> deltaPhiResiduals_BeforeClustering = new ArrayList<>();
    private static List<Double> deltaTimeResiduals_BeforeClustering = new ArrayList<>();

    // Lists for residuals after clustering
    private static List<Double> deltaZResiduals_Hits = new ArrayList<>();
    private static List<Double> deltaPhiResiduals_Hits = new ArrayList<>();
    private static List<Double> deltaTimeResiduals_Hits = new ArrayList<>();

    public ATOFHit_Reco_Cluster(String title) {
        super(title);
    }

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc") || !reader.getSchemaFactory().hasSchema("MC::True")) {
            System.err.println("Required schemas not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ATOFHit_Reco_Cluster demo = new ATOFHit_Reco_Cluster("Cluster Analysis Plots");
        demo.createPlots();

        System.out.println("\nCluster Size Summary:");
        for (Map.Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            System.out.printf("Clusters of size %d: %d\n", entry.getKey(), entry.getValue());
        }
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

            if (mcTrueBank.getRows() > 0) {
                // Calculate residuals before clustering
                System.out.println("Calculating residuals before clustering...");
                extractAndCompareTruth(mcTrueBank, hits, true);
            }

            // Separate bar and wedge hits
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

                processCluster(barLeft, barRight, wedgeHits, zBar, tBar);

                // After clustering, calculate residuals for clustered hits
                if (mcTrueBank.getRows() > 0) {
                    System.out.println("Calculating residuals after clustering...");
                    List<Hit> clusterHits = new ArrayList<>();
                    clusterHits.add(barLeft);
                    clusterHits.add(barRight);
                    clusterHits.addAll(wedgeHits);
                    extractAndCompareTruth(mcTrueBank, clusterHits, false);
                }
            }

            eventCount++;
        }
    }

    private static void processCluster(Hit barLeft, Hit barRight, List<Hit> wedgeHits, double zBar, double tBar) {
        List<Hit> clusterWedgeHits = new ArrayList<>();
        System.out.println("Clustering:");

        for (Hit wedgeHit : wedgeHits) {
            double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
            double deltaPhi = normalizePhi(barLeft.getPhi() - wedgeHit.getPhi());
            double deltaTime = Math.abs(tBar - wedgeHit.getTime());

            deltaZList.add(deltaZ);
            deltaPhiList.add(deltaPhi);
            deltaTimeList.add(deltaTime);

            System.out.printf("Wedge Hit: DeltaZ=%.2f, DeltaPhi=%.2f, DeltaTime=%.2f\n", deltaZ, deltaPhi, deltaTime);

            if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                clusterWedgeHits.add(wedgeHit);
            }
        }

        if (!clusterWedgeHits.isEmpty()) {
            int clusterSize = 2 + clusterWedgeHits.size();
            clusterSizes.add(clusterSize);
            clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

            double zCluster = calculateWeightedAverageZ(barLeft, barRight, clusterWedgeHits);
            double phiCluster = calculateWeightedAveragePhi(barLeft, barRight, clusterWedgeHits);
            double timeCluster = calculateWeightedAverageTime(barLeft, barRight, clusterWedgeHits);
            double energyCluster = calculateClusterEnergy(barLeft, barRight, clusterWedgeHits);

            zClusterList.add(zCluster);
            phiClusterList.add(phiCluster);
            timeClusterList.add(timeCluster);
            energyClusterList.add(energyCluster);

            System.out.printf("Cluster: ZCluster=%.2f, PhiCluster=%.2f, TimeCluster=%.2f, EnergyCluster=%.2f, Size=%d\n",
                    zCluster, phiCluster, timeCluster, energyCluster, clusterSize);
        }
    }

    private static void extractAndCompareTruth(Bank mcTrueBank, List<Hit> hits, boolean beforeClustering) {
        for (int i = 0; i < mcTrueBank.getRows(); i++) {
            double truthZ = mcTrueBank.getFloat("avgZ", i);
            double truthPhi = Math.atan2(mcTrueBank.getFloat("avgY", i), mcTrueBank.getFloat("avgX", i));
            double truthTime = mcTrueBank.getFloat("avgT", i);

            for (Hit hit : hits) {
                double residualZ = hit.zWedge() - truthZ;
                double residualPhi = normalizePhi(hit.getPhi() - truthPhi);
                double residualTime = hit.getTime() - truthTime;

                if (beforeClustering) {
                    deltaZResiduals_BeforeClustering.add(residualZ);
                    deltaPhiResiduals_BeforeClustering.add(residualPhi);
                    deltaTimeResiduals_BeforeClustering.add(residualTime);
                    System.out.printf("Before Clustering - Truth Z: %.2f, Reconstructed Z: %.2f, Residual Z: %.2f\n",
                            truthZ, hit.zWedge(), residualZ);
                } else {
                    deltaZResiduals_Hits.add(residualZ);
                    deltaPhiResiduals_Hits.add(residualPhi);
                    deltaTimeResiduals_Hits.add(residualTime);
                    System.out.printf("After Clustering - Truth Z: %.2f, Reconstructed Z: %.2f, Residual Z: %.2f\n",
                            truthZ, hit.zWedge(), residualZ);
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

    private static double calculateWeightedAverageZ(Hit barLeft, Hit barRight, List<Hit> clusterWedgeHits) {
        double weightedSumZ = barLeft.getAdc() * barLeft.zWedge() + barRight.getAdc() * barRight.zWedge();
        double totalWeight = barLeft.getAdc() + barRight.getAdc();

        for (Hit hit : clusterWedgeHits) {
            weightedSumZ += hit.getAdc() * hit.zWedge();
            totalWeight += hit.getAdc();
        }

        return totalWeight > 0 ? weightedSumZ / totalWeight : 0.0;
    }

    private static double calculateWeightedAveragePhi(Hit barLeft, Hit barRight, List<Hit> clusterWedgeHits) {
        double weightedSumPhi = barLeft.getAdc() * barLeft.getPhi() + barRight.getAdc() * barRight.getPhi();
        double totalWeight = barLeft.getAdc() + barRight.getAdc();

        for (Hit hit : clusterWedgeHits) {
            weightedSumPhi += hit.getAdc() * hit.getPhi();
            totalWeight += hit.getAdc();
        }

        return totalWeight > 0 ? weightedSumPhi / totalWeight : 0.0;
    }

    private static double calculateWeightedAverageTime(Hit barLeft, Hit barRight, List<Hit> clusterWedgeHits) {
        double weightedSumTime = barLeft.getAdc() * barLeft.getTime() + barRight.getAdc() * barRight.getTime();
        double totalWeight = barLeft.getAdc() + barRight.getAdc();

        for (Hit hit : clusterWedgeHits) {
            weightedSumTime += hit.getAdc() * hit.getTime();
            totalWeight += hit.getAdc();
        }

        return totalWeight > 0 ? weightedSumTime / totalWeight : 0.0;
    }

    private static double calculateClusterEnergy(Hit barLeft, Hit barRight, List<Hit> clusterWedgeHits) {
        double totalEnergy = barLeft.getAdc() + barRight.getAdc();

        for (Hit hit : clusterWedgeHits) {
            totalEnergy += hit.getAdc();
        }

        return totalEnergy;
    }

    private void createPlots() {
        System.out.printf("Delta Z Residuals Before Clustering (size: %d)\n", deltaZResiduals_BeforeClustering.size());
        System.out.printf("Delta Phi Residuals Before Clustering (size: %d)\n", deltaPhiResiduals_BeforeClustering.size());
        System.out.printf("Delta Time Residuals Before Clustering (size: %d)\n", deltaTimeResiduals_BeforeClustering.size());

        // Plot residuals before clustering
        createHistogramPlot("Delta Z (Truth vs Reconstructed Before Clustering)", "Delta Z (mm)", deltaZResiduals_BeforeClustering);
        createHistogramPlot("Delta Phi (Truth vs Reconstructed Before Clustering)", "Delta Phi (rad)", deltaPhiResiduals_BeforeClustering);
        createHistogramPlot("Delta Time (Truth vs Reconstructed Before Clustering)", "Delta Time (ns)", deltaTimeResiduals_BeforeClustering);

        // Plot residuals after clustering
        createHistogramPlot("Delta Z (Truth vs Reconstructed Hits)", "Delta Z (mm)", deltaZResiduals_Hits);
        createHistogramPlot("Delta Phi (Truth vs Reconstructed Hits)", "Delta Phi (rad)", deltaPhiResiduals_Hits);
        createHistogramPlot("Delta Time (Truth vs Reconstructed Hits)", "Delta Time (ns)", deltaTimeResiduals_Hits);

        // Plot clustering deltas
        createHistogramPlot("Delta Z (Clustering)", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Phi (Clustering)", "Delta Phi (rad)", deltaPhiList);
        createHistogramPlot("Delta Time (Clustering)", "Delta Time (ns)", deltaTimeList);
    }

    private void createHistogramPlot(String title, String xAxisLabel, List<Double> data) {
        if (data == null || data.isEmpty()) {
            System.out.printf("No data available for plot: %s\n", title);
            return;
        }

        HistogramDataset dataset = new HistogramDataset();
        double[] values = data.stream().mapToDouble(Double::doubleValue).toArray();
        dataset.addSeries(title, values, 50);
        JFreeChart histogram = ChartFactory.createHistogram(title, xAxisLabel, "Frequency", dataset, PlotOrientation.VERTICAL, true, true, false);
        ChartPanel chartPanel = new ChartPanel(histogram);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));
        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
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

*/




