//below this code previous version of code also works to produce residual related  plots 

package org.jlab.rec.atof.RESIDUALS_CALCULATIONS;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ATOFHit_Reco_Cluster extends JFrame {

    
    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 3.0; // mm
    private static final int N_WEDGE = 10;
    private static final double Z_THRESHOLD = 280.0;
    private static final double PHI_THRESHOLD = 0.1; // radians
    private static final double TIME_THRESHOLD = 1.7;

        private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

        private static List<Double> zClusterList = new ArrayList<>();
    private static List<Double> phiClusterList = new ArrayList<>();
    private static List<Double> timeClusterList = new ArrayList<>();
    private static List<Double> energyClusterList = new ArrayList<>();

        private static List<Double> deltaZResiduals = new ArrayList<>();
    private static List<Double> deltaPhiResiduals = new ArrayList<>();
    private static List<Double> deltaTimeResiduals = new ArrayList<>();

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
            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

            System.out.printf("\nProcessing Event #%d with %d hits.\n", eventCount, numHits);

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofAdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit);
                else if (hit.layer >= 10 && hit.layer <= 19) wedgeHits.add(hit);
            }

            if (barHits.size() >= 2) {
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);

                double zBar = VEFF * (barLeft.time - barRight.time) / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF,
                        barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);

                processCluster(barLeft, barRight, wedgeHits, zBar, tBar);
            }

            if (mcTrueBank.getRows() > 0) {
                extractAndCompareTruth(mcTrueBank, barHits);
            }
            eventCount++;
        }
    }

    private static void processCluster(Hit barLeft, Hit barRight, List<Hit> wedgeHits, double zBar, double tBar) {
        List<Hit> clusterWedgeHits = new ArrayList<>();
        System.out.println("Clustering:");

        for (Hit wedgeHit : wedgeHits) {
            double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
            double deltaPhi = normalizePhi(barLeft.phi - wedgeHit.phi);
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

    private static double calculateWeightedAverageZ(Hit barLeft, Hit barRight, List<Hit> clusterWedgeHits) {
        double weightedSumZ = barLeft.adc * barLeft.zWedge() + barRight.adc * barRight.zWedge();
        double totalWeight = barLeft.adc + barRight.adc;

        for (Hit hit : clusterWedgeHits) {
            weightedSumZ += hit.adc * hit.zWedge();
            totalWeight += hit.adc;
        }

        return totalWeight > 0 ? weightedSumZ / totalWeight : 0.0;
    }

    private static double calculateWeightedAveragePhi(Hit barLeft, Hit barRight, List<Hit> clusterWedgeHits) {
        double weightedSumPhi = barLeft.adc * barLeft.phi + barRight.adc * barRight.phi;
        double totalWeight = barLeft.adc + barRight.adc;

        for (Hit hit : clusterWedgeHits) {
            weightedSumPhi += hit.adc * hit.phi;
            totalWeight += hit.adc;
        }

        return totalWeight > 0 ? weightedSumPhi / totalWeight : 0.0;
    }

    private static double calculateWeightedAverageTime(Hit barLeft, Hit barRight, List<Hit> clusterWedgeHits) {
        double weightedSumTime = barLeft.adc * barLeft.time + barRight.adc * barRight.time;
        double totalWeight = barLeft.adc + barRight.adc;

        for (Hit hit : clusterWedgeHits) {
            weightedSumTime += hit.adc * hit.time;
            totalWeight += hit.adc;
        }

        return totalWeight > 0 ? weightedSumTime / totalWeight : 0.0;
    }

    private static double calculateClusterEnergy(Hit barLeft, Hit barRight, List<Hit> clusterWedgeHits) {
        double totalEnergy = barLeft.adc + barRight.adc;

        for (Hit hit : clusterWedgeHits) {
            totalEnergy += hit.adc;
        }

        return totalEnergy;
    }

    private static double normalizePhi(double phi) {
        while (phi > Math.PI) phi -= 2 * Math.PI;
        while (phi < -Math.PI) phi += 2 * Math.PI;
        return phi;
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

    private static void extractAndCompareTruth(Bank mcTrueBank, List<Hit> barHits) {
        for (int i = 0; i < mcTrueBank.getRows(); i++) {
            double truthZ = mcTrueBank.getFloat("avgZ", i);
            double truthPhi = Math.atan2(mcTrueBank.getFloat("avgY", i), mcTrueBank.getFloat("avgX", i));
            double truthTime = mcTrueBank.getFloat("avgT", i);

            for (Hit barHit : barHits) {
                double residualZ = barHit.zWedge() - truthZ;
                double residualPhi = normalizePhi(barHit.phi - truthPhi);
                double residualTime = barHit.time - truthTime;

                deltaZResiduals.add(residualZ);
                deltaPhiResiduals.add(residualPhi);
                deltaTimeResiduals.add(residualTime);

                System.out.printf("Truth Z: %.2f, Reconstructed Z: %.2f, Residual: %.2f\n", truthZ, barHit.zWedge(), residualZ);
                System.out.printf("Truth Phi: %.2f, Reconstructed Phi: %.2f, Residual Phi: %.2f\n", truthPhi, barHit.phi, residualPhi);
                System.out.printf("Truth Time: %.2f, Reconstructed Time: %.2f, Residual Time: %.2f\n", truthTime, barHit.time, residualTime);
            }
        }
    }

    private void createPlots() {
        createHistogramPlot("Delta Z (Truth vs Reconstructed)", "Delta Z (mm)", deltaZResiduals);
        createHistogramPlot("Delta Phi (Truth vs Reconstructed)", "Delta Phi (rad)", deltaPhiResiduals);
        createHistogramPlot("Delta Time (Truth vs Reconstructed)", "Delta Time (ns)", deltaTimeResiduals);

        createHistogramPlot("Delta Z (Clustering)", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Phi (Clustering)", "Delta Phi (rad)", deltaPhiList);
        createHistogramPlot("Delta Time (Clustering)", "Delta Time (ns)", deltaTimeList);

        createScatterPlot("Cluster Z vs Cluster Index", "Cluster Index", "Z (mm)", zClusterList);
        createScatterPlot("Cluster Phi vs Cluster Index", "Cluster Index", "Phi (rad)", phiClusterList);
        createScatterPlot("Cluster Time vs Cluster Index", "Cluster Index", "Time (ns)", timeClusterList);
        createScatterPlot("Cluster Energy vs Cluster Index", "Cluster Index", "Energy (ADC)", energyClusterList);
    }

    private void createHistogramPlot(String title, String xAxisLabel, List<Double> data) {
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

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<?> yData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < yData.size(); i++) series.add(i, ((Number) yData.get(i)).doubleValue());
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset, PlotOrientation.VERTICAL, true, true, false);
        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));
        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    static class Hit {
        int sector, layer, component, order, adc, ped;
        double time, phi;

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

        double zWedge() {
            int wedgeIndex = component % N_WEDGE;
            return (wedgeIndex - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;
        }
    }
}






/*
//plots residuals very well 

package org.jlab.rec.atof.RESIDUALS_CALCULATIONS.org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ATOFHit_Reco_Cluster extends JFrame {


    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10;
    private static final double Z_THRESHOLD = 280.0;
    private static final double PHI_THRESHOLD = 0.1; // radians
    private static final double TIME_THRESHOLD = 1.7;

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

    private static List<Double> zClusterList = new ArrayList<>();
    private static List<Double> phiClusterList = new ArrayList<>();
    private static List<Double> timeClusterList = new ArrayList<>();
    private static List<Double> energyClusterList = new ArrayList<>();

    private static List<Double> deltaZResiduals = new ArrayList<>();
    private static List<Double> deltaPhiResiduals = new ArrayList<>();
    private static List<Double> deltaTimeResiduals = new ArrayList<>();

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
            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

            System.out.printf("\nProcessing Event #%d with %d hits.\n", eventCount, numHits);

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofAdcBank, hitIndex);
                if (hit.layer == 0) {
                    barHits.add(hit);
                    System.out.printf("Bar Hit: Sector=%d, Layer=%d, Component=%d, Order=%d, ADC=%d, Time=%.2f, Ped=%d\n",
                            hit.sector, hit.layer, hit.component, hit.order, hit.adc, hit.time, hit.ped);
                } else if (hit.layer >= 10 && hit.layer <= 19) {
                    wedgeHits.add(hit);
                    System.out.printf("Wedge Hit: Sector=%d, Layer=%d, Component=%d, Order=%d, ADC=%d, Time=%.2f, Phi=%.2f, Z=%.2f, Ped=%d\n",
                            hit.sector, hit.layer, hit.component, hit.order, hit.adc, hit.time, hit.phi, hit.zWedge(), hit.ped);
                }
            }

            if (barHits.size() >= 2) {
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);

                double zBar = VEFF * (barLeft.time - barRight.time) / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF,
                        barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);

                System.out.printf("Bar-Level: ZBar=%.2f, TBar=%.2f\n", zBar, tBar);

                processCluster(barLeft, barRight, wedgeHits, zBar, tBar);
            }

            if (mcTrueBank.getRows() > 0) {
                extractAndCompareTruth(mcTrueBank, barHits);
            }
            eventCount++;
        }
    }

    private static void extractAndCompareTruth(Bank mcTrueBank, List<Hit> barHits) {
        for (int i = 0; i < mcTrueBank.getRows(); i++) {
            double truthZ = mcTrueBank.getFloat("avgZ", i);
            double truthPhi = Math.atan2(mcTrueBank.getFloat("avgY", i), mcTrueBank.getFloat("avgX", i));
            double truthTime = mcTrueBank.getFloat("avgT", i);

            for (Hit barHit : barHits) {
                double residualZ = barHit.zWedge() - truthZ;
                deltaZResiduals.add(residualZ);

                double residualPhi = normalizePhi(barHit.phi - truthPhi);
                deltaPhiResiduals.add(residualPhi);

                double residualTime = barHit.time - truthTime;
                deltaTimeResiduals.add(residualTime);

                System.out.printf("Truth Z: %.2f, Reconstructed Z: %.2f, Residual Z: %.2f\n", truthZ, barHit.zWedge(), residualZ);
                System.out.printf("Truth Phi: %.2f, Reconstructed Phi: %.2f, Residual Phi: %.2f\n", truthPhi, barHit.phi, residualPhi);
                System.out.printf("Truth Time: %.2f, Reconstructed Time: %.2f, Residual Time: %.2f\n", truthTime, barHit.time, residualTime);
            }
        }
    }

    private static void processCluster(Hit barLeft, Hit barRight, List<Hit> wedgeHits, double zBar, double tBar) {
        List<Hit> clusterWedgeHits = new ArrayList<>();
        System.out.println("Clustering:");

        for (Hit wedgeHit : wedgeHits) {
            double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
            double deltaPhi = normalizePhi(barLeft.phi - wedgeHit.phi);
            double deltaTime = Math.abs(tBar - wedgeHit.time);

            deltaZList.add(deltaZ);
            deltaPhiList.add(deltaPhi);
            deltaTimeList.add(deltaTime);

            if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                clusterWedgeHits.add(wedgeHit);
            }

            System.out.printf("Wedge Hit: ZWedge=%.2f, PhiWedge=%.2f, Time=%.2f, DeltaZ=%.2f, DeltaPhi=%.2f, DeltaTime=%.2f\n",
                    wedgeHit.zWedge(), wedgeHit.phi, wedgeHit.time, deltaZ, deltaPhi, deltaTime);
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

    private static double normalizePhi(double phi) {
        while (phi > Math.PI) phi -= 2 * Math.PI;
        while (phi < -Math.PI) phi += 2 * Math.PI;
        return phi;
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

    private void createPlots() {
        createHistogramPlot("Delta Z (Truth vs Reconstructed)", "Delta Z (mm)", deltaZResiduals);
        createHistogramPlot("Delta Phi (Truth vs Reconstructed)", "Delta Phi (rad)", deltaPhiResiduals);
        createHistogramPlot("Delta Time (Truth vs Reconstructed)", "Delta Time (ns)", deltaTimeResiduals);

        createHistogramPlot("Delta Z (Clustering)", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Phi (Clustering)", "Delta Phi (rad)", deltaPhiList);
        createHistogramPlot("Delta Time (Clustering)", "Delta Time (ns)", deltaTimeList);

        createScatterPlot("Cluster Z vs Cluster Index", "Cluster Index", "Z (mm)", zClusterList);
        createScatterPlot("Cluster Phi vs Cluster Index", "Cluster Index", "Phi (rad)", phiClusterList);
        createScatterPlot("Cluster Time vs Cluster Index", "Cluster Index", "Time (ns)", timeClusterList);
        createScatterPlot("Cluster Energy vs Cluster Index", "Cluster Index", "Energy (ADC)", energyClusterList);
    }

    private void createHistogramPlot(String title, String xAxisLabel, List<Double> data) {
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

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<?> yData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < yData.size(); i++) series.add(i, ((Number) yData.get(i)).doubleValue());
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset, PlotOrientation.VERTICAL, true, true, false);
        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));
        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    static class Hit {
        int sector, layer, component, order, adc, ped;
        double time, phi;

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

        double zWedge() {
            int wedgeIndex = component % N_WEDGE;
            return (wedgeIndex - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;
        }
    }
}


*/






