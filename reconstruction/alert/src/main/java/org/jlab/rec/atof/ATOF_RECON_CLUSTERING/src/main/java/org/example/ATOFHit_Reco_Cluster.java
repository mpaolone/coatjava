
package org.example;

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
    private static final double PHI_THRESHOLD = 0.3;
    private static final double TIME_THRESHOLD = 1.7;

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static List<Integer> barClusterSizes = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

    // New lists to store weighted average values
    private static List<Double> zClusterList = new ArrayList<>();
    private static List<Double> phiClusterList = new ArrayList<>();
    private static List<Double> timeClusterList = new ArrayList<>();
    private static List<Double> energyClusterList = new ArrayList<>();

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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Required schema not found in the HIPO file.");
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

            // Bar-only clustering: Requires exactly 2 hits
            if (barHits.size() == 2) {
                barClusterSizes.add(2);
            } else {
                barClusterSizes.add(null); // No valid bar cluster for this event
            }

            // Bar and Wedge Clustering
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

                    double zCluster = calculateWeightedAverageZ(barLeft, barRight, clusterWedgeHits);
                    double phiCluster = calculateWeightedAveragePhi(barLeft, barRight, clusterWedgeHits);
                    double timeCluster = calculateWeightedAverageTime(barLeft, barRight, clusterWedgeHits);
                    double energyCluster = calculateClusterEnergy(barLeft, barRight, clusterWedgeHits);

                    zClusterList.add(zCluster);
                    phiClusterList.add(phiCluster);
                    timeClusterList.add(timeCluster);
                    energyClusterList.add(energyCluster);

                    System.out.printf("Cluster Formed (Size: %d):\n", clusterSize);
                    System.out.printf("  Cluster Z: %.2f mm, Phi: %.2f rad, Time: %.2f ns, Total Energy: %.2f ADC\n", zCluster, phiCluster, timeCluster, energyCluster);
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

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);

        createClusterSizeOverlayPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes, barClusterSizes);
        createScatterPlot("Cluster Z vs Cluster Index", "Cluster Index", "Z (mm)", zClusterList);
        createScatterPlot("Cluster Phi vs Cluster Index", "Cluster Index", "Phi (rad)", phiClusterList);
        createScatterPlot("Cluster Time vs Cluster Index", "Cluster Index", "Time (ns)", timeClusterList);
        createScatterPlot("Cluster Energy vs Cluster Index", "Cluster Index", "Energy (ADC)", energyClusterList);
    }

    private void createHistogramPlot(String title, String xAxisLabel, List<Double> data) {
        HistogramDataset dataset = new HistogramDataset();
        double[] values = data.stream().mapToDouble(Double::doubleValue).toArray();
        dataset.addSeries(title, values, 50);
        JFreeChart histogram = ChartFactory.createHistogram(title, xAxisLabel, "Events", dataset, PlotOrientation.VERTICAL, true, true, false);
        ChartPanel chartPanel = new ChartPanel(histogram);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));
        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private void createClusterSizeOverlayPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> clusterSizes, List<Integer> barClusterSizes) {
        XYSeries clusterSeries = new XYSeries("Bar + Wedge Clusters");
        XYSeries barOnlySeries = new XYSeries("Bar Clusters Only");

        for (int i = 0; i < clusterSizes.size(); i++) {
            clusterSeries.add(i, clusterSizes.get(i));
            if (barClusterSizes.get(i) != null) { // Only add valid bar cluster sizes
                barOnlySeries.add(i, barClusterSizes.get(i));
            }
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(clusterSeries);
        dataset.addSeries(barOnlySeries);

        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset, PlotOrientation.VERTICAL, true, true, false);
        ChartPanel chartPanel = new ChartPanel(scatterPlot);
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
package org.example;
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
    private static final double PHI_THRESHOLD = 0.3;
    private static final double TIME_THRESHOLD = 1.7;

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

    // New lists to store weighted average values
    private static List<Double> zClusterList = new ArrayList<>();
    private static List<Double> phiClusterList = new ArrayList<>();
    private static List<Double> timeClusterList = new ArrayList<>();
    private static List<Double> energyClusterList = new ArrayList<>(); // Cluster energy

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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Required schema not found in the HIPO file.");
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

                    double zCluster = calculateWeightedAverageZ(barLeft, barRight, clusterWedgeHits);
                    double phiCluster = calculateWeightedAveragePhi(barLeft, barRight, clusterWedgeHits);
                    double timeCluster = calculateWeightedAverageTime(barLeft, barRight, clusterWedgeHits);
                    double energyCluster = calculateClusterEnergy(barLeft, barRight, clusterWedgeHits);

                    zClusterList.add(zCluster);
                    phiClusterList.add(phiCluster);
                    timeClusterList.add(timeCluster);
                    energyClusterList.add(energyCluster);

                    System.out.printf("Cluster Formed (Size: %d):\n", clusterSize);
                    System.out.printf("  Cluster Z: %.2f mm, Phi: %.2f rad, Time: %.2f ns, Total Energy: %.2f ADC\n", zCluster, phiCluster, timeCluster, energyCluster);

                    int wedgeCount = 1;
                    for (Hit wedgeHit : clusterWedgeHits) {
                        double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                        double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                        double deltaTime = Math.abs(tBar - wedgeHit.time);

                        printHitDetails(wedgeHit, "Wedge Hit #" + wedgeCount++);
                        System.out.printf("    Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", deltaZ, deltaPhi, deltaTime);
                    }
                    System.out.println("---- End of Cluster ----\n");
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

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
        createScatterPlot("Cluster Z vs Cluster Index", "Cluster Index", "Z (mm)", zClusterList);
        createScatterPlot("Cluster Phi vs Cluster Index", "Cluster Index", "Phi (rad)", phiClusterList);
        createScatterPlot("Cluster Time vs Cluster Index", "Cluster Index", "Time (ns)", timeClusterList);
        createScatterPlot("Cluster Energy vs Cluster Index", "Cluster Index", "Energy (ADC)", energyClusterList);
    }

    private void createHistogramPlot(String title, String xAxisLabel, List<Double> data) {
        HistogramDataset dataset = new HistogramDataset();
        double[] values = data.stream().mapToDouble(Double::doubleValue).toArray();
        dataset.addSeries(title, values, 50);
        JFreeChart histogram = ChartFactory.createHistogram(title, xAxisLabel, "Events", dataset, PlotOrientation.VERTICAL, true, true, false);
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

*/









/* does great jhob in plotting 

package org.example;
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
    private static final double Z_THRESHOLD = 60.0;
    private static final double PHI_THRESHOLD = 0.1;
    private static final double TIME_THRESHOLD = 1.3;

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

    // New lists to store weighted average values
    private static List<Double> zClusterList = new ArrayList<>();
    private static List<Double> phiClusterList = new ArrayList<>();
    private static List<Double> timeClusterList = new ArrayList<>();
    private static List<Double> energyClusterList = new ArrayList<>(); // Cluster energy

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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Required schema not found in the HIPO file.");
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

                    double zCluster = calculateWeightedAverageZ(barLeft, barRight, clusterWedgeHits);
                    double phiCluster = calculateWeightedAveragePhi(barLeft, barRight, clusterWedgeHits);
                    double timeCluster = calculateWeightedAverageTime(barLeft, barRight, clusterWedgeHits);
                    double energyCluster = calculateClusterEnergy(barLeft, barRight, clusterWedgeHits);

                    zClusterList.add(zCluster);
                    phiClusterList.add(phiCluster);
                    timeClusterList.add(timeCluster);
                    energyClusterList.add(energyCluster);

                    // Print cluster details
                    System.out.printf("Cluster #%d (Size: %d):\n", eventCount, clusterSize);
                    System.out.printf("  Weighted Z: %.2f mm\n", zCluster);
                    System.out.printf("  Weighted Phi: %.2f rad\n", phiCluster);
                    System.out.printf("  Weighted Time: %.2f ns\n", timeCluster);
                    System.out.printf("  Total Energy: %.2f ADC\n", energyCluster);

                    // Print delta values for each hit within the cluster
                    int wedgeCount = 1;
                    for (Hit wedgeHit : clusterWedgeHits) {
                        double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                        double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                        double deltaTime = Math.abs(tBar - wedgeHit.time);

                        System.out.printf("    Wedge Hit #%d -> Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", 
                                          wedgeCount++, deltaZ, deltaPhi, deltaTime);
                    }
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

    // Weighted average calculation functions
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

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
        createScatterPlot("Cluster Z vs Cluster Index", "Cluster Index", "Z (mm)", zClusterList);
        createScatterPlot("Cluster Phi vs Cluster Index", "Cluster Index", "Phi (rad)", phiClusterList);
        createScatterPlot("Cluster Time vs Cluster Index", "Cluster Index", "Time (ns)", timeClusterList);
        createScatterPlot("Cluster Energy vs Cluster Index", "Cluster Index", "Energy (ADC)", energyClusterList);
    }

    private void createHistogramPlot(String title, String xAxisLabel, List<Double> data) {
        HistogramDataset dataset = new HistogramDataset();
        double[] values = data.stream().mapToDouble(Double::doubleValue).toArray();
        dataset.addSeries(title, values, 50);
        JFreeChart histogram = ChartFactory.createHistogram(title, xAxisLabel, "Events", dataset, PlotOrientation.VERTICAL, true, true, false);
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

*/


















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
    private static final double Z_THRESHOLD = 60.0;
    private static final double PHI_THRESHOLD = 0.1;
    private static final double TIME_THRESHOLD = 1.3;

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

    // New lists to store weighted average values
    private static List<Double> zClusterList = new ArrayList<>();
    private static List<Double> phiClusterList = new ArrayList<>();
    private static List<Double> timeClusterList = new ArrayList<>();
    private static List<Double> energyClusterList = new ArrayList<>(); // Cluster energy

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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Required schema not found in the HIPO file.");
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
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

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

                    // Calculate weighted averages and energy for the cluster
                    zClusterList.add(calculateWeightedAverageZ(barLeft, barRight, clusterWedgeHits));
                    phiClusterList.add(calculateWeightedAveragePhi(barLeft, barRight, clusterWedgeHits));
                    timeClusterList.add(calculateWeightedAverageTime(barLeft, barRight, clusterWedgeHits));
                    energyClusterList.add(calculateClusterEnergy(barLeft, barRight, clusterWedgeHits));
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

    // Weighted average calculation functions
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

    private void createPlots() {
        createScatterPlot("Cluster Z vs Cluster Index", "Cluster Index", "Z (mm)", zClusterList);
        createScatterPlot("Cluster Phi vs Cluster Index", "Cluster Index", "Phi (rad)", phiClusterList);
        createScatterPlot("Cluster Time vs Cluster Index", "Cluster Index", "Time (ns)", timeClusterList);
        createScatterPlot("Cluster Energy vs Cluster Index", "Cluster Index", "Energy (ADC)", energyClusterList);
    }

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Double> yData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < yData.size(); i++) series.add(i, yData.get(i));
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

*/
    
















































































//older versions less relevant 


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
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN extends JFrame {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Number of distinct Z values for wedges
    private static final double Z_THRESHOLD = 150.0; // mm
    private static final double PHI_THRESHOLD = 0.1; // rad
    private static final double TIME_THRESHOLD = 1.0; // ns

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Cluster Analysis Plots");
        demo.createPlots();

        // Print summary of cluster sizes
        System.out.println("\nCluster Size Summary:");
        for (Map.Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            System.out.printf("Clusters of size %d: %d\n", entry.getKey(), entry.getValue());
        }
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofAdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit); // Bar hits
                else if (hit.layer >= 10 && hit.layer <= 19) wedgeHits.add(hit); // Wedge hits
            }

            // Process only if there are at least 2 bar hits
            if (barHits.size() >= 2) {
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);

                // Calculate ZBar and TBar for the bar
                double zBar = VEFF * (barLeft.time - barRight.time) / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF, barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);

                // Cross-match each bar with all possible wedges and find valid wedges
                List<Hit> clusterWedgeHits = new ArrayList<>();
                for (Hit wedgeHit : wedgeHits) {
                    // Ensure the wedge hits have matching Z and Phi values
                    double zWedge = wedgeHit.zWedge();
                    double deltaZ = Math.abs(zBar - zWedge);
                    double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        deltaZList.add(deltaZ);
                        deltaPhiList.add(deltaPhi);
                        deltaTimeList.add(deltaTime);
                        clusterWedgeHits.add(wedgeHit);
                    }
                }

                // Only form a cluster if at least one wedge hit is within the thresholds
                if (!clusterWedgeHits.isEmpty()) {
                    int clusterSize = 2 + clusterWedgeHits.size();
                    clusterSizes.add(clusterSize);
                    clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

                    // Print detailed cluster information
                    System.out.printf("Cluster Formed (Size %d):\n", clusterSize);
                    System.out.printf("  ZBar: %.2f mm\n", zBar);
                    System.out.printf("  TBar: %.2f ns\n", tBar);
                    System.out.printf("  Bar Hit #1 -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                            barLeft.sector, barLeft.layer, barLeft.component, barLeft.order, barLeft.adc, barLeft.time, barLeft.phi);
                    System.out.printf("  Bar Hit #2 -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                            barRight.sector, barRight.layer, barRight.component, barRight.order, barRight.adc, barRight.time, barRight.phi);

                    System.out.println("  Wedge Hits:");
                    int wedgeCount = 1;
                    for (Hit wedgeHit : clusterWedgeHits) {
                        double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                        double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                        double deltaTime = Math.abs(tBar - wedgeHit.time);

                        System.out.printf("    Wedge Hit #%d -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, ZWedge: %.2f mm, Phi: %.2f rad)\n",
                                wedgeCount++, wedgeHit.sector, wedgeHit.layer, wedgeHit.component, wedgeHit.order, wedgeHit.adc, wedgeHit.time, wedgeHit.zWedge(), wedgeHit.phi);
                        System.out.printf("      Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", deltaZ, deltaPhi, deltaTime);
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

    private void createPlots() {
        if (!deltaZList.isEmpty()) {
            createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        } else {
            System.out.println("No data available for Delta Z Distribution.");
        }

        if (!deltaTimeList.isEmpty()) {
            createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        } else {
            System.out.println("No data available for Delta Time Distribution.");
        }

        if (!deltaPhiList.isEmpty()) {
            createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        } else {
            System.out.println("No data available for Delta Phi Distribution.");
        }

        if (!clusterSizes.isEmpty()) {
            createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
        } else {
            System.out.println("No data available for Cluster Size vs Event Index.");
        }
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

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < xData.size(); i++) series.add(i, xData.get(i));
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

*/







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
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN extends JFrame {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Number of distinct Z values for wedges
    private static final double Z_THRESHOLD = 150.0; // mm
    private static final double PHI_THRESHOLD = 0.1; // rad
    private static final double TIME_THRESHOLD = 1.0; // ns

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Cluster Analysis Plots");
        demo.createPlots();

        // Print summary of cluster sizes
        System.out.println("\nCluster Size Summary:");
        for (Map.Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            System.out.printf("Clusters of size %d: %d\n", entry.getKey(), entry.getValue());
        }
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofAdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit); // Bar hits
                else if (hit.layer >= 10 && hit.layer <= 19) wedgeHits.add(hit); // Wedge hits
            }

            if (barHits.size() >= 2) {
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);

                // Calculate ZBar, TBar, and PhiBar for the bar
                double zBar = VEFF * (barLeft.time - barRight.time) / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF, barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);
                double phiBar = barLeft.phi;

                // Enforce matching Z and Phi for wedges
                List<Hit> clusterWedgeHits = new ArrayList<>();
                for (Hit wedgeHit : wedgeHits) {
                    double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                    double deltaPhi = Math.abs(phiBar - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        // Check if the Z and Phi of this wedge hit are consistent with the cluster
                        if (clusterWedgeHits.isEmpty() || (clusterWedgeHits.get(0).zWedge() == wedgeHit.zWedge() && clusterWedgeHits.get(0).phi == wedgeHit.phi)) {
                            clusterWedgeHits.add(wedgeHit);
                        }
                    }
                }

                if (clusterWedgeHits.size() > 0) {
                    int clusterSize = 2 + clusterWedgeHits.size();
                    clusterSizes.add(clusterSize);
                    clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

                    // Print detailed cluster information
                    System.out.printf("Cluster Formed (Size %d):\n", clusterSize);
                    System.out.printf("  ZBar: %.2f mm\n", zBar);
                    System.out.printf("  TBar (min of TBarLeft and TBarRight): %.2f ns\n", tBar);
                    System.out.printf("  Bar Hit #1 -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                            barLeft.sector, barLeft.layer, barLeft.component, barLeft.order, barLeft.adc, barLeft.time, barLeft.phi);
                    System.out.printf("  Bar Hit #2 -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                            barRight.sector, barRight.layer, barRight.component, barRight.order, barRight.adc, barRight.time, barRight.phi);

                    System.out.println("  Wedge Hits:");
                    int wedgeCount = 1;
                    for (Hit wedgeHit : clusterWedgeHits) {
                        double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                        double deltaPhi = Math.abs(phiBar - wedgeHit.phi);
                        double deltaTime = Math.abs(tBar - wedgeHit.time);

                        System.out.printf("    Wedge Hit #%d -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, ZWedge: %.2f mm, Phi: %.2f rad)\n",
                                wedgeCount++, wedgeHit.sector, wedgeHit.layer, wedgeHit.component, wedgeHit.order, wedgeHit.adc, wedgeHit.time, wedgeHit.zWedge(), wedgeHit.phi);
                        System.out.printf("      Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", deltaZ, deltaPhi, deltaTime);
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

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
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

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < xData.size(); i++) series.add(i, xData.get(i));
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


*/




/* near correct 
package org.example;

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

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN extends JFrame {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Number of distinct Z values for wedges
    private static final double Z_THRESHOLD = 280.0; // mm
    private static final double PHI_THRESHOLD = 0.5; // rad
    private static final double TIME_THRESHOLD = 1.7; // ns

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Cluster Analysis Plots");
        demo.createPlots();

        // Print summary of cluster sizes
        System.out.println("\nCluster Size Summary:");
        for (Map.Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            System.out.printf("Clusters of size %d: %d\n", entry.getKey(), entry.getValue());
        }
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofAdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit); // Bar hits
                else if (hit.layer >= 10 && hit.layer <= 19) wedgeHits.add(hit); // Wedge hits
            }

            // Process only if there are at least 2 bar hits
            if (barHits.size() >= 2) {
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);

                // Calculate ZBar for the bar
                double zBar = VEFF * (barLeft.time - barRight.time) / 2;

                // Calculate ZLeft and ZRight
                double zLeft = zBar - BAR_LENGTH / 2;
                double zRight = zBar + BAR_LENGTH / 2;

                // Calculate TBarLeft and TBarRight based on ZLeft and ZRight
                double tBarLeft = barLeft.time - zLeft / VEFF;
                double tBarRight = barRight.time - zRight / VEFF;

                // Calculate TBar as the minimum of TBarLeft and TBarRight
                double tBar = Math.min(tBarLeft, tBarRight);

                // Cross-match each bar with all possible wedges and find valid wedges
                List<Hit> clusterWedgeHits = new ArrayList<>();
                for (Hit wedgeHit : wedgeHits) {
                    // Calculate delta values relative to the bar position and timing
                    double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                    double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    // Add delta values to respective lists for plotting
                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    // Add wedge hit to cluster if it meets the thresholds
                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        clusterWedgeHits.add(wedgeHit);
                    }
                }

                // Only form a cluster if at least one wedge hit is within the thresholds
                if (clusterWedgeHits.size() > 0) {
                    int clusterSize = 2 + clusterWedgeHits.size();
                    clusterSizes.add(clusterSize);
                    clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

                    // Print detailed cluster information
                    System.out.printf("Cluster Formed (Size %d):\n", clusterSize);
                    System.out.printf("  ZBar: %.2f mm\n", zBar);
                    System.out.printf("  TBar (min of TBarLeft and TBarRight): %.2f ns\n", tBar);
                    System.out.printf("  Bar Hit #1 -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                            barLeft.sector, barLeft.layer, barLeft.component, barLeft.order, barLeft.adc, barLeft.time, barLeft.phi);
                    System.out.printf("  Bar Hit #2 -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                            barRight.sector, barRight.layer, barRight.component, barRight.order, barRight.adc, barRight.time, barRight.phi);

                    System.out.println("  Wedge Hits:");
                    int wedgeCount = 1;
                    for (Hit wedgeHit : clusterWedgeHits) {
                        double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                        double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                        double deltaTime = Math.abs(tBar - wedgeHit.time);

                        System.out.printf("    Wedge Hit #%d -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, ZWedge: %.2f mm, Phi: %.2f rad)\n",
                                wedgeCount++, wedgeHit.sector, wedgeHit.layer, wedgeHit.component, wedgeHit.order, wedgeHit.adc, wedgeHit.time, wedgeHit.zWedge(), wedgeHit.phi);
                        System.out.printf("      Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", deltaZ, deltaPhi, deltaTime);
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

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
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

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < xData.size(); i++) series.add(i, xData.get(i));
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

*/





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
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN  extends JFrame {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Number of distinct Z values for wedges
    private static final double Z_THRESHOLD = 150.0; // mm
    private static final double PHI_THRESHOLD = 0.1; // rad
    private static final double TIME_THRESHOLD = 1.0; // ns

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Cluster Analysis Plots");
        demo.createPlots();

        // Print summary of cluster sizes
        System.out.println("\nCluster Size Summary:");
        for (Map.Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            System.out.printf("Clusters of size %d: %d\n", entry.getKey(), entry.getValue());
        }
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofAdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit); // Bar hits
                else if (hit.layer >= 10 && hit.layer <= 19) wedgeHits.add(hit); // Wedge hits
            }

            if (barHits.size() >= 2) {
                // Calculate ZBar and TBar based on two bar hits
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);
                double deltaT = barLeft.time - barRight.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF, barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);

                // Match bar with wedge hits; only form clusters if wedge hits meet conditions
                List<Hit> clusterWedgeHits = new ArrayList<>();

                for (Hit wedgeHit : wedgeHits) {
                    double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                    double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    // Store delta values for histogram plotting
                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    // If wedge hit meets threshold criteria, include it in the cluster
                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        clusterWedgeHits.add(wedgeHit);
                    }
                }

                // Only print the cluster if it includes at least one wedge hit
                if (!clusterWedgeHits.isEmpty()) {
                    int clusterSize = 2 + clusterWedgeHits.size(); // Two bar hits plus wedge hits
                    clusterSizes.add(clusterSize);
                    clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

                    // Print cluster details
                    System.out.printf("Cluster Formed (Size %d):\n", clusterSize);
                    System.out.printf("  Bar Hit #1 -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, ZBar: %.2f mm, Phi: %.2f rad)\n",
                            barLeft.sector, barLeft.layer, barLeft.component, barLeft.order, barLeft.adc, barLeft.time, zBar, barLeft.phi);
                    System.out.printf("  Bar Hit #2 -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, ZBar: %.2f mm, Phi: %.2f rad)\n",
                            barRight.sector, barRight.layer, barRight.component, barRight.order, barRight.adc, barRight.time, zBar, barRight.phi);
                    System.out.println("  Wedge Hits:");
                    int wedgeIndex = 1;
                    for (Hit wedgeHit : clusterWedgeHits) {
                        System.out.printf("    Wedge Hit #%d -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, ZWedge: %.2f mm, Phi: %.2f rad)\n",
                                wedgeIndex++, wedgeHit.sector, wedgeHit.layer, wedgeHit.component, wedgeHit.order, wedgeHit.adc, wedgeHit.time, wedgeHit.zWedge(), wedgeHit.phi);
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

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
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

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < xData.size(); i++) series.add(i, xData.get(i));
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

*/












/* //works good 
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.statistics.HistogramDataset;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN  extends JFrame {

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

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processHitsGlobally(reader);
        reader.close();

      ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Global Cluster Analysis Plots");
        demo.createPlots();

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

                    System.out.printf("Cluster Formed (Size %d):\n", clusterSize);
                    System.out.printf("  Bar Hits:\n");
                    System.out.printf("    Bar Hit 1 -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                                      barHit1.sector, barHit1.layer, barHit1.component, barHit1.order, barHit1.adc, barHit1.time, barHit1.phi);
                    System.out.printf("    Bar Hit 2 -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                                      barHit2.sector, barHit2.layer, barHit2.component, barHit2.order, barHit2.adc, barHit2.time, barHit2.phi);
                    System.out.printf("  Combined Bar Z and T: ZBar: %.2f mm, TBar: %.2f ns\n", zBar, tBar);

                    System.out.println("  Wedge Hits:");
                    int wedgeCount = 1;
                    for (Hit wedgeHit : clusterWedgeHits) {
                        double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                        double deltaPhi = Math.abs(barHit1.phi - wedgeHit.phi);
                        double deltaTime = Math.abs(tBar - wedgeHit.time);
                        System.out.printf("    Wedge Hit %d -> (Sector: %d, Layer: %d, Component: %d, ZWedge: %.2f mm, Phi: %.2f rad, Time: %.2f ns), Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n",
                                          wedgeCount++, wedgeHit.sector, wedgeHit.layer, wedgeHit.component, wedgeHit.zWedge(), wedgeHit.phi, wedgeHit.time,
                                          deltaZ, deltaPhi, deltaTime);
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

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
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

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < xData.size(); i++) series.add(i, xData.get(i));
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
*/








/*
package org.example;
// Updated code to enforce bar hits on the same bar

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.statistics.HistogramDataset;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN  extends JFrame {

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

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processHitsGlobally(reader);
        reader.close();

       ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Global Cluster Analysis Plots");
        demo.createPlots();

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

        // Process clusters across bar and wedge hits, ensuring bar hits come from the same bar
        for (int i = 0; i < barHits.size(); i++) {
            Hit barHit1 = barHits.get(i);

            for (int j = i + 1; j < barHits.size(); j++) {
                Hit barHit2 = barHits.get(j);

                // Only consider bar hits from the same bar
                if (barHit1.sector == barHit2.sector &&
                    barHit1.layer == barHit2.layer &&
                    barHit1.component == barHit2.component) {

                    // Calculate combined Z and T for the bar
                    double deltaT = barHit1.time - barHit2.time;
                    double zBar = VEFF * deltaT / 2;
                    double tBar = Math.min(barHit1.time - (zBar - BAR_LENGTH / 2) / VEFF, 
                                           barHit2.time - (zBar + BAR_LENGTH / 2) / VEFF);

                    // Initialize cluster with the two bar hits
                    int clusterSize = 2;

                    // Check proximity to each wedge hit
                    for (Hit wedgeHit : wedgeHits) {
                        double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                        double deltaPhi = Math.abs(barHit1.phi - wedgeHit.phi);
                        double deltaTime = Math.abs(tBar - wedgeHit.time);

                        // Add values to plotting lists
                        deltaZList.add(deltaZ);
                        deltaPhiList.add(deltaPhi);
                        deltaTimeList.add(deltaTime);

                        // Form clusters if criteria met
                        if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                            clusterSize++;
                            clusterSizes.add(clusterSize);
                            clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

                            // Print details for clusters
                            System.out.printf("Cluster Formed (Size %d):\n", clusterSize);
                            System.out.printf("  Bar Hits:\n");
                            System.out.printf("    Bar Hit 1 -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                                              barHit1.sector, barHit1.layer, barHit1.component, barHit1.order, barHit1.adc, barHit1.time, barHit1.phi);
                            System.out.printf("    Bar Hit 2 -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                                              barHit2.sector, barHit2.layer, barHit2.component, barHit2.order, barHit2.adc, barHit2.time, barHit2.phi);
                            System.out.printf("  Combined Bar Z and T: ZBar: %.2f mm, TBar: %.2f ns\n", zBar, tBar);
                            System.out.printf("  Wedge Hit -> (Sector: %d, Layer: %d, Component: %d, ZWedge: %.2f mm, Phi: %.2f rad, Time: %.2f ns), Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n",
                                              wedgeHit.sector, wedgeHit.layer, wedgeHit.component, wedgeHit.zWedge(), wedgeHit.phi, wedgeHit.time, deltaZ, deltaPhi, deltaTime);
                        }
                    }
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

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
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

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < xData.size(); i++) series.add(i, xData.get(i));
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
*/





/*

package org.example;


// Updated code to enforce bar hits on the same bar

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.statistics.HistogramDataset;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN  extends JFrame {

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

    public  ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processHitsGlobally(reader);
        reader.close();

       ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN  demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Global Cluster Analysis Plots");
        demo.createPlots();

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

        // Process clusters across bar and wedge hits, ensuring bar hits come from the same bar
        for (int i = 0; i < barHits.size(); i++) {
            Hit barHit1 = barHits.get(i);

            for (int j = i + 1; j < barHits.size(); j++) {
                Hit barHit2 = barHits.get(j);

                // Only consider bar hits from the same bar
                if (barHit1.sector == barHit2.sector &&
                    barHit1.layer == barHit2.layer &&
                    barHit1.component == barHit2.component) {

                    // Calculate combined Z and T for the bar
                    double deltaT = barHit1.time - barHit2.time;
                    double zBar = VEFF * deltaT / 2;
                    double tBar = Math.min(barHit1.time - (zBar - BAR_LENGTH / 2) / VEFF, 
                                           barHit2.time - (zBar + BAR_LENGTH / 2) / VEFF);

                    // Initialize cluster with the two bar hits
                    int clusterSize = 2;

                    // Check proximity to each wedge hit
                    for (Hit wedgeHit : wedgeHits) {
                        double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                        double deltaPhi = Math.abs(barHit1.phi - wedgeHit.phi);
                        double deltaTime = Math.abs(tBar - wedgeHit.time);

                        // Add values to plotting lists
                        deltaZList.add(deltaZ);
                        deltaPhiList.add(deltaPhi);
                        deltaTimeList.add(deltaTime);

                        // Form clusters if criteria met
                        if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                            clusterSize++;
                            clusterSizes.add(clusterSize);
                            clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

                            // Print details for clusters
                            System.out.printf("Cluster Formed (Size %d):\n", clusterSize);
                            System.out.printf("  Bar Hits:\n");
                            System.out.printf("    Bar Hit 1 -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                                              barHit1.sector, barHit1.layer, barHit1.component, barHit1.order, barHit1.adc, barHit1.time, barHit1.phi);
                            System.out.printf("    Bar Hit 2 -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                                              barHit2.sector, barHit2.layer, barHit2.component, barHit2.order, barHit2.adc, barHit2.time, barHit2.phi);
                            System.out.printf("  Combined Bar Z and T: ZBar: %.2f mm, TBar: %.2f ns\n", zBar, tBar);
                            System.out.printf("  Wedge Hit -> (Sector: %d, Layer: %d, Component: %d, ZWedge: %.2f mm, Phi: %.2f rad, Time: %.2f ns), Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n",
                                              wedgeHit.sector, wedgeHit.layer, wedgeHit.component, wedgeHit.zWedge(), wedgeHit.phi, wedgeHit.time, deltaZ, deltaPhi, deltaTime);
                        }
                    }
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

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
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

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < xData.size(); i++) series.add(i, xData.get(i));
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


*/




/*
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.statistics.HistogramDataset;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN  extends JFrame {

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
    private static boolean detailedPrinted = false;

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processHitsGlobally(reader);
        reader.close();

       ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN  demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Global Cluster Analysis Plots");
        demo.createPlots();

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

        // Process clusters across all bar and wedge hits
        for (int i = 0; i < barHits.size(); i++) {
            for (int j = i + 1; j < barHits.size(); j++) {
                Hit barLeft = barHits.get(i);
                Hit barRight = barHits.get(j);

                // Calculate ZBar and TBar based on two bar hits
                double deltaT = barLeft.time - barRight.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF, barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);

                // Match each wedge across all bars
                int clusterSize = 2;
                List<Hit> matchedWedges = new ArrayList<>();
                for (Hit wedgeHit : wedgeHits) {
                    double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                    double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    // Add values to full range plots
                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    // Form clusters if criteria met
                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        clusterSize++;
                        matchedWedges.add(wedgeHit);
                    }
                }
                
                clusterSizes.add(clusterSize);
                clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

                // Print only detailed information for a large cluster, such as size >= 207
                if (!detailedPrinted && clusterSize >= 207) {
                    System.out.printf("\nCluster Formed (Size %d):\n", clusterSize);
                    System.out.println("  Bar Hits:");
                    System.out.printf("    Bar Hit -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                            barLeft.sector, barLeft.layer, barLeft.component, barLeft.order, barLeft.adc, barLeft.time, barLeft.phi);
                    System.out.printf("    Bar Hit -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                            barRight.sector, barRight.layer, barRight.component, barRight.order, barRight.adc, barRight.time, barRight.phi);
                    System.out.printf("  Combined Bar Z and T: ZBar: %.2f mm, TBar: %.2f ns\n", zBar, tBar);

                    System.out.println("  Wedge Hits:");
                    for (Hit wedgeHit : matchedWedges) {
                        double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                        double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                        double deltaTime = Math.abs(tBar - wedgeHit.time);
                        System.out.printf("    Wedge Hit -> (Sector: %d, Layer: %d, Component: %d, ZWedge: %.2f mm, Phi: %.2f rad, Time: %.2f ns)\n",
                                wedgeHit.sector, wedgeHit.layer, wedgeHit.component, wedgeHit.zWedge(), wedgeHit.phi, wedgeHit.time);
                        System.out.printf("      Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", deltaZ, deltaPhi, deltaTime);
                    }
                    System.out.println("---- End of Cluster ----\n");
                    detailedPrinted = true;
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

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
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

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < xData.size(); i++) series.add(i, xData.get(i));
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

*/







/* broad clustering all bars and all wedges
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.statistics.HistogramDataset;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN  extends JFrame {

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

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processHitsGlobally(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN  demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Global Cluster Analysis Plots");
        demo.createPlots();

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

        // Process clusters across all bar and wedge hits
        for (int i = 0; i < barHits.size(); i++) {
            for (int j = i + 1; j < barHits.size(); j++) {
                Hit barLeft = barHits.get(i);
                Hit barRight = barHits.get(j);

                // Calculate ZBar and TBar based on two bar hits
                double deltaT = barLeft.time - barRight.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF, barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);

                // Match each wedge across all bars
                int clusterSize = 2;
                for (Hit wedgeHit : wedgeHits) {
                    double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                    double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    // Add values to full range plots
                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    // Form clusters if criteria met
                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        clusterSize++;
                        clusterSizes.add(clusterSize);
                        clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

                        System.out.printf("Cluster Formed (Size %d) -> Bar Hits: ZBar: %.2f mm, TBar: %.2f ns, Wedge Hit: (Sector: %d, Layer: %d, Component: %d, ZWedge: %.2f mm, Phi: %.2f rad, Time: %.2f ns), Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n",
                                clusterSize, zBar, tBar, wedgeHit.sector, wedgeHit.layer, wedgeHit.component, wedgeHit.zWedge(), wedgeHit.phi, wedgeHit.time, deltaZ, deltaPhi, deltaTime);
                    }
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

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
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

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < xData.size(); i++) series.add(i, xData.get(i));
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


*/





/*
//works based on per event logic best code 
package org.example;

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

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN extends JFrame {

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

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Cluster Analysis Plots");
        demo.createPlots();

        // Print summary of cluster sizes
        System.out.println("\nCluster Size Summary:");
        for (Map.Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            System.out.printf("Clusters of size %d: %d\n", entry.getKey(), entry.getValue());
        }
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofAdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit); // Bar hits
                else if (hit.layer >= 10 && hit.layer <= 19) wedgeHits.add(hit); // Wedge hits
            }

            if (barHits.size() >= 2) {
                // Calculate ZBar and TBar based on two bar hits
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);
                double deltaT = barLeft.time - barRight.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF, barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);

                // Incrementally build clusters with additional wedge hits
                List<Hit> clusterWedges = new ArrayList<>();
                for (int i = 0; i < wedgeHits.size(); i++) {
                    Hit wedgeHit = wedgeHits.get(i);

                    // Calculate deltas between bar and each wedge hit
                    double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                    double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    // Add values to full range plots
                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    // Form clusters with incremental wedge hits if criteria met
                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        clusterWedges.add(wedgeHit);
                        int clusterSize = 2 + clusterWedges.size(); // 2 bar hits + wedge count
                        clusterSizes.add(clusterSize);
                        clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

                        System.out.printf("Cluster Formed (Size %d)\n", clusterSize);
                        System.out.println("  Bar Hits:");
                        System.out.printf("    (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                                barLeft.sector, barLeft.layer, barLeft.component, barLeft.order, barLeft.adc, barLeft.time, barLeft.phi);
                        System.out.printf("    (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                                barRight.sector, barRight.layer, barRight.component, barRight.order, barRight.adc, barRight.time, barRight.phi);

                        System.out.printf("  Combined Bar Info -> ZBar: %.2f mm, TBar: %.2f ns\n", zBar, tBar);

                        int wedgeCount = 1;
                        for (Hit clusteredWedge : clusterWedges) {
                            System.out.printf("  Wedge Hit #%d -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, ZWedge: %.2f mm, Phi: %.2f rad)\n",
                                    wedgeCount, clusteredWedge.sector, clusteredWedge.layer, clusteredWedge.component, clusteredWedge.order, clusteredWedge.adc, clusteredWedge.time, clusteredWedge.zWedge(), clusteredWedge.phi);
                            System.out.printf("    Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n",
                                    Math.abs(zBar - clusteredWedge.zWedge()), Math.abs(barLeft.phi - clusteredWedge.phi), Math.abs(tBar - clusteredWedge.time));
                            wedgeCount++;
                        }
                    }
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

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
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

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < xData.size(); i++) series.add(i, xData.get(i));
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
package org.example;

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

public class package org.example;

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

public class AdvancedClusterAnalysis extends JFrame {

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

    public AdvancedClusterAnalysis(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        AdvancedClusterAnalysis demo = new AdvancedClusterAnalysis("Cluster Analysis Plots");
        demo.createPlots();

        // Print summary of cluster sizes
        System.out.println("\nCluster Size Summary:");
        for (Map.Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            System.out.printf("Clusters of size %d: %d\n", entry.getKey(), entry.getValue());
        }
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofAdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit); // Bar hits
                else if (hit.layer >= 10 && hit.layer <= 19) wedgeHits.add(hit); // Wedge hits
            }

            if (barHits.size() >= 2) {
                // Calculate ZBar and TBar based on two bar hits
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);
                double deltaT = barLeft.time - barRight.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF, barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);

                // Incrementally build clusters with additional wedge hits
                List<Hit> clusterWedges = new ArrayList<>();
                for (int i = 0; i < wedgeHits.size(); i++) {
                    Hit wedgeHit = wedgeHits.get(i);

                    // Calculate deltas between bar and each wedge hit
                    double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                    double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    // Add values to full range plots
                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    // Form clusters with incremental wedge hits if criteria met
                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        clusterWedges.add(wedgeHit);
                        int clusterSize = 2 + clusterWedges.size(); // 2 bar hits + wedge count
                        clusterSizes.add(clusterSize);
                        clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

                        System.out.printf("Cluster Formed (Size %d)\n", clusterSize);
                        System.out.println("  Bar Hits:");
                        System.out.printf("    (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                                barLeft.sector, barLeft.layer, barLeft.component, barLeft.order, barLeft.adc, barLeft.time, barLeft.phi);
                        System.out.printf("    (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                                barRight.sector, barRight.layer, barRight.component, barRight.order, barRight.adc, barRight.time, barRight.phi);

                        System.out.printf("  Combined Bar Info -> ZBar: %.2f mm, TBar: %.2f ns\n", zBar, tBar);

                        int wedgeCount = 1;
                        for (Hit clusteredWedge : clusterWedges) {
                            System.out.printf("  Wedge Hit #%d -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, ZWedge: %.2f mm, Phi: %.2f rad)\n",
                                    wedgeCount, clusteredWedge.sector, clusteredWedge.layer, clusteredWedge.component, clusteredWedge.order, clusteredWedge.adc, clusteredWedge.time, clusteredWedge.zWedge(), clusteredWedge.phi);
                            System.out.printf("    Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n",
                                    Math.abs(zBar - clusteredWedge.zWedge()), Math.abs(barLeft.phi - clusteredWedge.phi), Math.abs(tBar - clusteredWedge.time));
                            wedgeCount++;
                        }
                    }
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

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
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

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < xData.size(); i++) series.add(i, xData.get(i));
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
 extends JFrame {

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

    public AdvancedClusterAnalysis(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        AdvancedClusterAnalysis demo = new AdvancedClusterAnalysis("Cluster Analysis Plots");
        demo.createPlots();

        // Print summary of cluster sizes
        System.out.println("\nCluster Size Summary:");
        for (Map.Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            System.out.printf("Clusters of size %d: %d\n", entry.getKey(), entry.getValue());
        }
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofAdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit); // Bar hits
                else if (hit.layer >= 10 && hit.layer <= 19) wedgeHits.add(hit); // Wedge hits
            }

            if (barHits.size() >= 2) {
                // Calculate ZBar and TBar based on two bar hits
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);
                double deltaT = barLeft.time - barRight.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF, barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);

                // Incrementally build clusters with additional wedge hits
                List<Hit> clusterWedges = new ArrayList<>();
                for (int i = 0; i < wedgeHits.size(); i++) {
                    Hit wedgeHit = wedgeHits.get(i);

                    // Calculate deltas between bar and each wedge hit
                    double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                    double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    // Add values to full range plots
                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    // Form clusters with incremental wedge hits if criteria met
                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        clusterWedges.add(wedgeHit);
                        int clusterSize = 2 + clusterWedges.size(); // 2 bar hits + wedge count
                        clusterSizes.add(clusterSize);
                        clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

                        System.out.printf("Cluster Formed (Size %d)\n", clusterSize);
                        System.out.println("  Bar Hits:");
                        System.out.printf("    (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                                barLeft.sector, barLeft.layer, barLeft.component, barLeft.order, barLeft.adc, barLeft.time, barLeft.phi);
                        System.out.printf("    (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad)\n",
                                barRight.sector, barRight.layer, barRight.component, barRight.order, barRight.adc, barRight.time, barRight.phi);

                        System.out.printf("  Combined Bar Info -> ZBar: %.2f mm, TBar: %.2f ns\n", zBar, tBar);

                        int wedgeCount = 1;
                        for (Hit clusteredWedge : clusterWedges) {
                            System.out.printf("  Wedge Hit #%d -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, ZWedge: %.2f mm, Phi: %.2f rad)\n",
                                    wedgeCount, clusteredWedge.sector, clusteredWedge.layer, clusteredWedge.component, clusteredWedge.order, clusteredWedge.adc, clusteredWedge.time, clusteredWedge.zWedge(), clusteredWedge.phi);
                            System.out.printf("    Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n",
                                    Math.abs(zBar - clusteredWedge.zWedge()), Math.abs(barLeft.phi - clusteredWedge.phi), Math.abs(tBar - clusteredWedge.time));
                            wedgeCount++;
                        }
                    }
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

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
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

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < xData.size(); i++) series.add(i, xData.get(i));
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

*/






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
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN extends JFrame {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Number of distinct Z values for wedges
    private static final double Z_THRESHOLD = 150.0; // mm
    private static final double PHI_THRESHOLD = 0.1; // rad
    private static final double TIME_THRESHOLD = 1.0; // ns

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

       ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN  demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Cluster Analysis Plots");
        demo.createPlots();

        // Print summary of cluster sizes
        System.out.println("\nCluster Size Summary:");
        for (Map.Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            System.out.printf("Clusters of size %d: %d\n", entry.getKey(), entry.getValue());
        }
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofAdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit); // Bar hits
                else if (hit.layer >= 10 && hit.layer <= 19) wedgeHits.add(hit); // Wedge hits
            }

            if (barHits.size() >= 2) {
                // Calculate ZBar and TBar based on two bar hits
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);
                double deltaT = barLeft.time - barRight.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF, barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);

                // Cross-match each bar with all possible wedges
                for (int i = 0; i < wedgeHits.size(); i++) {
                    Hit wedgeHit = wedgeHits.get(i);

                    // Calculate deltas between bar and each wedge hit
                    double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                    double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    // Add values to full range plots
                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    // Form clusters if criteria met
                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        int clusterSize = 2 + (i + 1); // 2 bar hits + wedge count
                        clusterSizes.add(clusterSize);
                        clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

                        System.out.printf("Cluster Formed (Size %d)\n", clusterSize);
                        System.out.printf("  Bar Hit #1 -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, ZBar: %.2f mm, Phi: %.2f rad)\n",
                                barLeft.sector, barLeft.layer, barLeft.component, barLeft.order, barLeft.adc, barLeft.time, zBar, barLeft.phi);
                        System.out.printf("  Bar Hit #2 -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, ZBar: %.2f mm, Phi: %.2f rad)\n",
                                barRight.sector, barRight.layer, barRight.component, barRight.order, barRight.adc, barRight.time, zBar, barRight.phi);
                        System.out.printf("  Wedge Hit #%d -> (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, ZWedge: %.2f mm, Phi: %.2f rad)\n",
                                i + 1, wedgeHit.sector, wedgeHit.layer, wedgeHit.component, wedgeHit.order, wedgeHit.adc, wedgeHit.time, wedgeHit.zWedge(), wedgeHit.phi);
                        System.out.printf("    Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", deltaZ, deltaPhi, deltaTime);
                    }
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

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
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

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < xData.size(); i++) series.add(i, xData.get(i));
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

*/





/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN  extends JFrame {

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

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc") || !reader.getSchemaFactory().hasSchema("REC::Particle")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        // Print cluster summary
        System.out.println("\nCluster Size Summary:");
        clusterSizeCounts.forEach((size, count) -> System.out.printf("Clusters of size %d: %d\n", size, count));
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Bank particleBank = new Bank(reader.getSchemaFactory().getSchema("REC::Particle"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);
            event.read(particleBank);

            int numHits = atofAdcBank.getRows();
            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofAdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit); // Bar hits
                else if (hit.layer >= 10 && hit.layer <= 19) wedgeHits.add(hit); // Wedge hits
            }

            if (barHits.size() >= 2) {
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);
                double deltaT = barLeft.time - barRight.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF, barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);

                int wedgeIndex = 1;
                for (Hit wedgeHit : wedgeHits) {
                    double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                    double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        int clusterSize = 2 + wedgeIndex;
                        clusterSizes.add(clusterSize);
                        clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

                        System.out.printf("Cluster Formed (Size %d) -> Bar Hits: (ZBar: %.2f mm, TBar: %.2f ns), Wedge Hit #%d: (Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, ZWedge: %.2f mm), Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n",
                                clusterSize, zBar, tBar, wedgeIndex, wedgeHit.sector, wedgeHit.layer, wedgeHit.component, wedgeHit.order, wedgeHit.adc, wedgeHit.time, wedgeHit.zWedge(), deltaZ, deltaPhi, deltaTime);
                        wedgeIndex++;
                    }
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

*/







/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.statistics.HistogramDataset;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN  extends JFrame {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Number of wedges per bar
    private static final double Z_THRESHOLD = 280.0; // mm
    private static final double PHI_THRESHOLD = 0.2; // rad
    private static final double TIME_THRESHOLD = 2.0; // ns

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static List<Double> protonAngles = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc") || !reader.getSchemaFactory().hasSchema("REC::Particle")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN  demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Cluster Analysis Plots");
        demo.createPlots();

        // Print summary of cluster sizes
        System.out.println("\nCluster Size Summary:");
        for (Map.Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            System.out.printf("Clusters of size %d: %d\n", entry.getKey(), entry.getValue());
        }
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Bank particleBank = new Bank(reader.getSchemaFactory().getSchema("REC::Particle"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);
            event.read(particleBank);

            double protonAngle = getProtonAngle(particleBank);
            if (protonAngle != -1) protonAngles.add(protonAngle);

            int numHits = atofAdcBank.getRows();
            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofAdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit); // Bar hits
                else if (hit.layer >= 10 && hit.layer <= 19) wedgeHits.add(hit); // Wedge hits
            }

            if (barHits.size() >= 2) {
                // Calculate ZBar and TBar based on two bar hits
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);
                double deltaT = barLeft.time - barRight.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF, barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);

                // Incrementally build clusters with additional wedge hits
                for (int i = 0; i < wedgeHits.size(); i++) {
                    Hit wedgeHit = wedgeHits.get(i);

                    // Calculate deltas between bar and each wedge hit
                    double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                    double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    // Add values to full range plots
                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    // Form clusters with incremental wedge hits if criteria met
                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        int clusterSize = 2 + (i + 1); // 2 bar hits + wedge count
                        clusterSizes.add(clusterSize);
                        clusterSizeCounts.put(clusterSize, clusterSizeCounts.getOrDefault(clusterSize, 0) + 1);

                        System.out.printf("Cluster Formed (Size %d) -> Bar Hits: (ZBar: %.2f mm, TBar: %.2f ns), Wedge Hit #%d: (ZWedge: %.2f mm, Time: %.2f ns), Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n",
                                clusterSize, zBar, tBar, i + 1, wedgeHit.zWedge(), wedgeHit.time, deltaZ, deltaPhi, deltaTime);
                    }
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

    private static double getProtonAngle(Bank particleBank) {
        for (int i = 0; i < particleBank.getRows(); i++) {
            int pid = particleBank.getInt("pid", i);
            if (pid == 2212) { // Proton
                float px = particleBank.getFloat("px", i);
                float py = particleBank.getFloat("py", i);
                float pz = particleBank.getFloat("pz", i);
                return Math.toDegrees(Math.atan2(Math.sqrt(px * px + py * py), pz));
            }
        }
        return -1;
    }

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
        createScatterPlotForAngleVsDeltaZ("Delta Z vs Proton Angle", "Proton Angle (degrees)", "Delta Z (mm)", protonAngles, deltaZList);
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

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < xData.size(); i++) series.add(i, xData.get(i));
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset, PlotOrientation.VERTICAL, true, true, false);
        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));
        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private void createScatterPlotForAngleVsDeltaZ(String title, String xAxisLabel, String yAxisLabel, List<Double> angles, List<Double> deltaZ) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < Math.min(angles.size(), deltaZ.size()); i++) series.add(angles.get(i), deltaZ.get(i));
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
*/







/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.statistics.HistogramDataset;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN extends JFrame {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm
    private static final double PHI_THRESHOLD = 0.1; // rad
    private static final double TIME_THRESHOLD = 1.0; // ns

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static List<Double> protonAngles = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc") || !reader.getSchemaFactory().hasSchema("REC::Particle")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Cluster Analysis Plots");
        demo.createPlots();

        // Print summary of cluster sizes
        System.out.println("\nCluster Size Summary:");
        for (Map.Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            System.out.printf("Clusters of size %d: %d\n", entry.getKey(), entry.getValue());
        }
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Bank particleBank = new Bank(reader.getSchemaFactory().getSchema("REC::Particle"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);
            event.read(particleBank);

            double protonAngle = getProtonAngle(particleBank);
            if (protonAngle != -1) protonAngles.add(protonAngle);

            int numHits = atofAdcBank.getRows();
            int clustersInEvent = 0;
            List<Hit> barHits = new ArrayList<>();
            List<Hit> validWedgeHits = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofAdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit);
                else if (hit.layer >= 10 && hit.layer <= 19) validWedgeHits.add(hit);
            }

            if (barHits.size() >= 2) {
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);

                double deltaT = barLeft.time - barRight.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF, barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);

                for (Hit wedgeHit : validWedgeHits) {
                    double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                    double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    // Add to full range delta plots
                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    // Check if this forms a valid cluster with the current bars
                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        clustersInEvent++;
                        System.out.printf("Cluster Formed -> Bar Hits: (ZBar: %.2f mm, TBar: %.2f ns), Wedge Hit: (ZWedge: %.2f mm, Time: %.2f ns) Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n",
                                zBar, tBar, wedgeHit.zWedge(), wedgeHit.time, deltaZ, deltaPhi, deltaTime);
                    }
                }
            }

            if (clustersInEvent >= 3) {
                clusterSizes.add(clustersInEvent);
                clusterSizeCounts.put(clustersInEvent, clusterSizeCounts.getOrDefault(clustersInEvent, 0) + 1);
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

    private static double getProtonAngle(Bank particleBank) {
        for (int i = 0; i < particleBank.getRows(); i++) {
            int pid = particleBank.getInt("pid", i);
            if (pid == 2212) { // Proton
                float px = particleBank.getFloat("px", i);
                float py = particleBank.getFloat("py", i);
                float pz = particleBank.getFloat("pz", i);
                return Math.toDegrees(Math.atan2(Math.sqrt(px * px + py * py), pz));
            }
        }
        return -1;
    }

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlotForClusterSize("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
        createScatterPlotForAngleVsDeltaZ("Delta Z vs Proton Angle", "Proton Angle (degrees)", "Delta Z (mm)", protonAngles, deltaZList);
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

    private void createScatterPlotForClusterSize(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < xData.size(); i++) {
            series.add(i, xData.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset, PlotOrientation.VERTICAL, true, true, false);
        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));
        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private void createScatterPlotForAngleVsDeltaZ(String title, String xAxisLabel, String yAxisLabel, List<Double> angles, List<Double> deltaZ) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < Math.min(angles.size(), deltaZ.size()); i++) {
            series.add(angles.get(i), deltaZ.get(i));
        }

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


*/







/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.statistics.HistogramDataset;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN extends JFrame {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm
    private static final double PHI_THRESHOLD = 0.1; // rad
    private static final double TIME_THRESHOLD = 1.0; // ns

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static List<Double> protonAngles = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc") || !reader.getSchemaFactory().hasSchema("REC::Particle")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Cluster Analysis Plots");
        demo.createPlots();

        // Print summary of cluster sizes
        System.out.println("\nCluster Size Summary:");
        for (Map.Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            System.out.printf("Clusters of size %d: %d\n", entry.getKey(), entry.getValue());
        }
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Bank particleBank = new Bank(reader.getSchemaFactory().getSchema("REC::Particle"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);
            event.read(particleBank);

            double protonAngle = getProtonAngle(particleBank);
            if (protonAngle != -1) protonAngles.add(protonAngle);

            int numHits = atofAdcBank.getRows();
            int clustersInEvent = 0;
            List<Hit> barHits = new ArrayList<>();
            List<Hit> validWedgeHits = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofAdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit);
                else if (hit.layer >= 10 && hit.layer <= 19) validWedgeHits.add(hit);
            }

            if (barHits.size() >= 2) {
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);

                double deltaT = barLeft.time - barRight.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF, barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);

                for (Hit wedgeHit : validWedgeHits) {
                    double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                    double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        clustersInEvent++;
                    }
                }
            }

            if (clustersInEvent >= 3) {
                clusterSizes.add(clustersInEvent);
                clusterSizeCounts.put(clustersInEvent, clusterSizeCounts.getOrDefault(clustersInEvent, 0) + 1);
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

    private static double getProtonAngle(Bank particleBank) {
        for (int i = 0; i < particleBank.getRows(); i++) {
            int pid = particleBank.getInt("pid", i);
            if (pid == 2212) { // Proton
                float px = particleBank.getFloat("px", i);
                float py = particleBank.getFloat("py", i);
                float pz = particleBank.getFloat("pz", i);
                return Math.toDegrees(Math.atan2(Math.sqrt(px * px + py * py), pz));
            }
        }
        return -1;
    }

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlotForClusterSize("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
        createScatterPlotForAngleVsDeltaZ("Delta Z vs Proton Angle", "Proton Angle (degrees)", "Delta Z (mm)", protonAngles, deltaZList);
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

    private void createScatterPlotForClusterSize(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < xData.size(); i++) {
            series.add(i, xData.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset, PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));
        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private void createScatterPlotForAngleVsDeltaZ(String title, String xAxisLabel, String yAxisLabel, List<Double> xData, List<Double> yData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < Math.min(xData.size(), yData.size()); i++) {
            series.add(xData.get(i), yData.get(i));
        }

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
*/






/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.statistics.HistogramDataset;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN extends JFrame {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm
    private static final double PHI_THRESHOLD = 0.1; // rad
    private static final double TIME_THRESHOLD = 1.0; // ns

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static List<Double> protonAngles = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc") || !reader.getSchemaFactory().hasSchema("REC::Particle")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Cluster Analysis Plots");
        demo.createPlots();

        // Print summary of cluster sizes
        System.out.println("\nCluster Size Summary:");
        for (Map.Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            System.out.printf("Clusters of size %d: %d\n", entry.getKey(), entry.getValue());
        }
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Bank particleBank = new Bank(reader.getSchemaFactory().getSchema("REC::Particle"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);
            event.read(particleBank);

            double protonAngle = getProtonAngle(particleBank);
            if (protonAngle != -1) protonAngles.add(protonAngle);

            int numHits = atofAdcBank.getRows();
            int clustersInEvent = 0;
            List<Hit> barHits = new ArrayList<>();
            List<Hit> validWedgeHits = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofAdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit);
                else if (hit.layer >= 10 && hit.layer <= 19) validWedgeHits.add(hit);
            }

            if (barHits.size() >= 2) {
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);

                double deltaT = barLeft.time - barRight.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF, barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);

                for (Hit wedgeHit : validWedgeHits) {
                    double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                    double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        clustersInEvent++;
                    }
                }
            }

            if (clustersInEvent >= 3) {
                clusterSizes.add(clustersInEvent);
                clusterSizeCounts.put(clustersInEvent, clusterSizeCounts.getOrDefault(clustersInEvent, 0) + 1);
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

    private static double getProtonAngle(Bank particleBank) {
        for (int i = 0; i < particleBank.getRows(); i++) {
            int pid = particleBank.getInt("pid", i);
            if (pid == 2212) { // Proton
                float px = particleBank.getFloat("px", i);
                float py = particleBank.getFloat("py", i);
                float pz = particleBank.getFloat("pz", i);
                return Math.toDegrees(Math.atan2(Math.sqrt(px * px + py * py), pz));
            }
        }
        return -1;
    }

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
        createScatterPlot("Delta Z vs Proton Angle", "Proton Angle (degrees)", "Delta Z (mm)", protonAngles, deltaZList);
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

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData, List<Double> yData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < Math.min(xData.size(), yData.size()); i++) {
            series.add(xData.get(i), yData.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset, PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));
        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Double> xData, List<Double> yData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < Math.min(xData.size(), yData.size()); i++) {
            series.add(xData.get(i), yData.get(i));
        }

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
        int sector;
        int layer;
        int component;
        int order;
        int adc;
        double time;
        double phi;

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
            return (component % N_WEDGE - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;
        }
    }
}

*/











/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.statistics.HistogramDataset;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN extends JFrame {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm
    private static final double PHI_THRESHOLD = 0.1; // rad
    private static final double TIME_THRESHOLD = 1.0; // ns

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static List<Double> protonAngles = new ArrayList<>();
    private static Map<Integer, Integer> clusterSizeCounts = new HashMap<>();

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc") || !reader.getSchemaFactory().hasSchema("REC::Particle")) {
            System.err.println("Required schema not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Cluster Analysis Plots");
        demo.createPlots();

        // Print summary of cluster sizes
        System.out.println("\nCluster Size Summary:");
        for (Map.Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            System.out.printf("Clusters of size %d: %d\n", entry.getKey(), entry.getValue());
        }
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Bank particleBank = new Bank(reader.getSchemaFactory().getSchema("REC::Particle"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);
            event.read(particleBank);

            double protonAngle = getProtonAngle(particleBank);
            if (protonAngle != -1) protonAngles.add(protonAngle);

            int numHits = atofAdcBank.getRows();
            int clustersInEvent = 0;
            List<Hit> barHits = new ArrayList<>();
            List<Hit> validWedgeHits = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                Hit hit = createHit(atofAdcBank, hitIndex);
                if (hit.layer == 0) barHits.add(hit);
                else if (hit.layer >= 10 && hit.layer <= 19) validWedgeHits.add(hit);
            }

            if (barHits.size() >= 2) {
                Hit barLeft = barHits.get(0);
                Hit barRight = barHits.get(1);

                double deltaT = barLeft.time - barRight.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeft.time - (zBar - BAR_LENGTH / 2) / VEFF, barRight.time - (zBar + BAR_LENGTH / 2) / VEFF);

                for (Hit wedgeHit : validWedgeHits) {
                    double deltaZ = Math.abs(zBar - wedgeHit.zWedge());
                    double deltaPhi = Math.abs(barLeft.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        clustersInEvent++;
                    }
                }
            }

            if (clustersInEvent >= 3) {
                clusterSizes.add(clustersInEvent);
                clusterSizeCounts.put(clustersInEvent, clusterSizeCounts.getOrDefault(clustersInEvent, 0) + 1);
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

    private static double getProtonAngle(Bank particleBank) {
        for (int i = 0; i < particleBank.getRows(); i++) {
            int pid = particleBank.getInt("pid", i);
            if (pid == 2212) { // Proton
                float px = particleBank.getFloat("px", i);
                float py = particleBank.getFloat("py", i);
                float pz = particleBank.getFloat("pz", i);
                return Math.toDegrees(Math.atan2(Math.sqrt(px * px + py * py), pz));
            }
        }
        return -1;
    }

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
        createScatterPlot("Delta Z vs Proton Angle", "Proton Angle (degrees)", "Delta Z (mm)", protonAngles, deltaZList);
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

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData, List<Double> yData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < Math.min(xData.size(), yData.size()); i++) {
            series.add(xData.get(i), yData.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset, PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));
        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Double> xData, List<Double> yData) {
        XYSeries series = new XYSeries("Data");
        for (int i = 0; i < Math.min(xData.size(), yData.size()); i++) {
            series.add(xData.get(i), yData.get(i));
        }

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
        int sector;
        int layer;
        int component;
        int order;
        int adc;
        double time;
        double phi;

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
            return (component % N_WEDGE - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;
        }
    }
}


*/








/* 
//it wont plot delta Z vs Proton ANgle/

package org.example;

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
import java.util.List;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN extends JFrame {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 200.0; // mm, threshold for clustering
    private static final double PHI_THRESHOLD = 0.2; // rad, threshold for clustering
    private static final double TIME_THRESHOLD = 2.0; // ns, threshold for clustering

    private static List<Double> deltaZFullRangeList = new ArrayList<>();
    private static List<Double> deltaPhiFullRangeList = new ArrayList<>();
    private static List<Double> deltaTimeFullRangeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static List<Double> protonAngles = new ArrayList<>();
    private static List<Double> deltaZProtonAngleList = new ArrayList<>();

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
        super(title);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc") || !reader.getSchemaFactory().hasSchema("REC::Particle")) {
            System.err.println("Schema ATOF::adc or REC::Particle not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Cluster Analysis Plots");
        demo.createPlots();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Bank recParticleBank = new Bank(reader.getSchemaFactory().getSchema("REC::Particle"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);
            event.read(recParticleBank);

            int numHits = atofAdcBank.getRows();
            int clustersInEvent = 0;
            double protonAngle = 0.0;

            // Calculate proton angle from REC::Particle bank (assuming proton PID is 2212)
            for (int i = 0; i < recParticleBank.getRows(); i++) {
                int pid = recParticleBank.getInt("pid", i);
                if (pid == 2212) { // Proton PID
                    double px = recParticleBank.getFloat("px", i);
                    double py = recParticleBank.getFloat("py", i);
                    double pz = recParticleBank.getFloat("pz", i);
                    protonAngle = Math.toDegrees(Math.acos(pz / Math.sqrt(px * px + py * py + pz * pz)));
                    protonAngles.add(protonAngle);
                    break;
                }
            }

            Hit barLeftHit = null;
            Hit barRightHit = null;

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int sector = atofAdcBank.getByte("sector", hitIndex);
                int layer = atofAdcBank.getByte("layer", hitIndex);
                int component = atofAdcBank.getShort("component", hitIndex);
                int order = atofAdcBank.getByte("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                float time = atofAdcBank.getFloat("time", hitIndex);
                int pedestal = atofAdcBank.getShort("ped", hitIndex);
                double phi = 2 * Math.PI * component / 60;

                System.out.printf("Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                        sector, layer, component, order, adc, time, pedestal, phi);

                if (layer == 0) {
                    if (order == 0) {
                        barLeftHit = new Hit(sector, layer, component, order, adc, time, phi, pedestal);
                    } else {
                        barRightHit = new Hit(sector, layer, component, order, adc, time, phi, pedestal);
                    }

                    if (barLeftHit != null && barRightHit != null && barLeftHit.getComponent() == barRightHit.getComponent()) {
                        double deltaT = barLeftHit.getTime() - barRightHit.getTime();
                        double zBar = VEFF * deltaT / 2;
                        double tBar = Math.min(barLeftHit.getTime() - (zBar - BAR_LENGTH / 2) / VEFF,
                                barRightHit.getTime() - (zBar + BAR_LENGTH / 2) / VEFF);

                        barLeftHit = null;
                        barRightHit = null;

                        for (int wedgeIndex = 0; wedgeIndex < numHits; wedgeIndex++) {
                            int wedgeLayer = atofAdcBank.getByte("layer", wedgeIndex);
                            if (wedgeLayer >= 10 && wedgeLayer <= 19) {
                                int wedgeComponent = atofAdcBank.getShort("component", wedgeIndex);
                                double wedgeTime = atofAdcBank.getFloat("time", wedgeIndex);
                                double wedgePhi = phi;
                                int wedgeIndexWithinBar = wedgeComponent % N_WEDGE;
                                double zWedge = (wedgeIndexWithinBar - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;

                                double deltaZ = Math.abs(zBar - zWedge);
                                double deltaPhi = Math.abs(phi - wedgePhi);
                                double deltaTime = Math.abs(tBar - wedgeTime);

                                // Record deltas for full range plots
                                deltaZFullRangeList.add(deltaZ);
                                deltaPhiFullRangeList.add(deltaPhi);
                                deltaTimeFullRangeList.add(deltaTime);

                                // Apply thresholds for valid clusters
                                if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                                    clustersInEvent++;
                                    deltaZProtonAngleList.add(deltaZ);
                                }
                            }
                        }
                    }
                }
            }
            // Only clustered events are added to the cluster size plot
            if (clustersInEvent >= 3) {
                clusterSizes.add(clustersInEvent);
            }
        }
    }

    private void createPlots() {
        createHistogramPlot("Delta Z Full Range", "Delta Z (mm)", deltaZFullRangeList);
        createHistogramPlot("Delta Time Full Range", "Delta Time (ns)", deltaTimeFullRangeList);
        createHistogramPlot("Delta Phi Full Range", "Delta Phi (rad)", deltaPhiFullRangeList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
        createAnglePlot("Delta Z vs Proton Angle", "Proton Angle (degrees)", "Delta Z (mm)", protonAngles, deltaZProtonAngleList);
    }

    private void createHistogramPlot(String title, String xAxisLabel, List<Double> data) {
        HistogramDataset dataset = new HistogramDataset();
        double[] values = data.stream().mapToDouble(Double::doubleValue).toArray();
        dataset.addSeries(title, values, 50);

        JFreeChart histogram = ChartFactory.createHistogram(title, xAxisLabel, "Frequency", dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(histogram);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> data) {
        XYSeries series = new XYSeries("Cluster Size");
        for (int i = 0; i < data.size(); i++) {
            series.add(i, data.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private void createAnglePlot(String title, String xAxisLabel, String yAxisLabel, List<Double> xData, List<Double> yData) {
        XYSeries series = new XYSeries("Delta Z vs Angle");
        for (int i = 0; i < xData.size(); i++) {
            series.add(xData.get(i), yData.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart anglePlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(anglePlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    static class Hit {
        private int sector;
        private int layer;
        private int component;
        private int order;
        private int adc;
        private double time;
        private double phi;
        private int pedestal;

        public Hit(int sector, int layer, int component, int order, int adc, double time, double phi, int pedestal) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.time = time;
            this.phi = phi;
            this.pedestal = pedestal;
        }

        public int getComponent() { return component; }
        public double getTime() { return time; }
    }
}
*/









/*
//it extracts all hits with sector layer components and order for both bar and wedges

package org.example;

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
import java.util.List;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN extends JFrame {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm, threshold for clustering
    private static final double PHI_THRESHOLD = 0.1; // rad, threshold for clustering
    private static final double TIME_THRESHOLD = 1.0; // ns, threshold for clustering

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static List<Double> protonAngles = new ArrayList<>();
    private static List<Double> deltaZProtonAngleList = new ArrayList<>();

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
        super(title);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc") || !reader.getSchemaFactory().hasSchema("REC::Particle")) {
            System.err.println("Schema ATOF::adc or REC::Particle not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Cluster Analysis Plots");
        demo.createPlots();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Bank recParticleBank = new Bank(reader.getSchemaFactory().getSchema("REC::Particle"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);
            event.read(recParticleBank);

            int numHits = atofAdcBank.getRows();
            int clustersInEvent = 0;

            Hit barLeftHit = null;
            Hit barRightHit = null;
            double protonAngle = 0.0;

            // Calculate proton angle from REC::Particle bank (assuming proton PID is 2212)
            for (int i = 0; i < recParticleBank.getRows(); i++) {
                int pid = recParticleBank.getInt("pid", i);
                if (pid == 2212) { // Proton PID
                    double px = recParticleBank.getFloat("px", i);
                    double py = recParticleBank.getFloat("py", i);
                    double pz = recParticleBank.getFloat("pz", i);
                    protonAngle = Math.toDegrees(Math.acos(pz / Math.sqrt(px * px + py * py + pz * pz)));
                    protonAngles.add(protonAngle);
                    break;
                }
            }

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int sector = atofAdcBank.getByte("sector", hitIndex);
                int layer = atofAdcBank.getByte("layer", hitIndex);
                int component = atofAdcBank.getShort("component", hitIndex);
                int order = atofAdcBank.getByte("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                float time = atofAdcBank.getFloat("time", hitIndex);
                int pedestal = atofAdcBank.getShort("ped", hitIndex);
                double phi = 2 * Math.PI * component / 60;

                System.out.printf("Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                        sector, layer, component, order, adc, time, pedestal, phi);

                if (layer == 0) {
                    if (order == 0) {
                        barLeftHit = new Hit(sector, layer, component, order, adc, time, phi, pedestal);
                    } else {
                        barRightHit = new Hit(sector, layer, component, order, adc, time, phi, pedestal);
                    }

                    if (barLeftHit != null && barRightHit != null && barLeftHit.getComponent() == barRightHit.getComponent()) {
                        double deltaT = barLeftHit.getTime() - barRightHit.getTime();
                        double zBar = VEFF * deltaT / 2;
                        double tBar = Math.min(barLeftHit.getTime() - (zBar - BAR_LENGTH / 2) / VEFF,
                                barRightHit.getTime() - (zBar + BAR_LENGTH / 2) / VEFF);

                        barLeftHit = null;
                        barRightHit = null;

                        for (int wedgeIndex = 0; wedgeIndex < numHits; wedgeIndex++) {
                            int wedgeLayer = atofAdcBank.getByte("layer", wedgeIndex);
                            if (wedgeLayer >= 10 && wedgeLayer <= 19) {
                                int wedgeComponent = atofAdcBank.getShort("component", wedgeIndex);
                                double wedgeTime = atofAdcBank.getFloat("time", wedgeIndex);
                                double wedgePhi = phi;
                                int wedgeIndexWithinBar = wedgeComponent % N_WEDGE;
                                double zWedge = (wedgeIndexWithinBar - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;

                                double deltaZ = Math.abs(zBar - zWedge);
                                double deltaPhi = Math.abs(phi - wedgePhi);
                                double deltaTime = Math.abs(tBar - wedgeTime);

                                if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                                    clustersInEvent++;
                                    deltaZList.add(deltaZ);
                                    deltaPhiList.add(deltaPhi);
                                    deltaTimeList.add(deltaTime);
                                    deltaZProtonAngleList.add(deltaZ);
                                }
                            }
                        }
                    }
                }
            }
            clusterSizes.add(clustersInEvent);
        }
    }

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes);
        createAnglePlot("Delta Z vs Proton Angle", "Proton Angle (degrees)", "Delta Z (mm)", protonAngles, deltaZProtonAngleList);
    }

    private void createHistogramPlot(String title, String xAxisLabel, List<Double> data) {
        HistogramDataset dataset = new HistogramDataset();
        double[] values = data.stream().mapToDouble(Double::doubleValue).toArray();
        dataset.addSeries(title, values, 50);

        JFreeChart histogram = ChartFactory.createHistogram(title, xAxisLabel, "Frequency", dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(histogram);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> data) {
        XYSeries series = new XYSeries("Cluster Size");
        for (int i = 0; i < data.size(); i++) {
            series.add(i, data.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private void createAnglePlot(String title, String xAxisLabel, String yAxisLabel, List<Double> xData, List<Double> yData) {
        XYSeries series = new XYSeries("Delta Z vs Angle");
        for (int i = 0; i < xData.size(); i++) {
            series.add(xData.get(i), yData.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart anglePlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(anglePlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    static class Hit {
        private int sector;
        private int layer;
        private int component;
        private int order;
        private int adc;
        private double time;
        private double phi;
        private int pedestal;

        public Hit(int sector, int layer, int component, int order, int adc, double time, double phi, int pedestal) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.time = time;
            this.phi = phi;
            this.pedestal = pedestal;
        }

        public int getComponent() { return component; }
        public double getTime() { return time; }
    }
}
*/






/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm, threshold for delta Z
    private static final double PHI_THRESHOLD = 0.1; // rad, threshold for delta Phi
    private static final double TIME_THRESHOLD = 1.0; // ns, threshold for delta Time

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Double> clusterSizes = new ArrayList<>(); // Changed to List<Double> for compatibility
    private static List<Double> eventIndices = new ArrayList<>(); // Changed to List<Double> for compatibility
    private static List<Double> angleList = new ArrayList<>(); // List to store angles for plotting

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN  demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN();
        demo.createPlots();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();
        int eventIndex = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            System.out.println("\nProcessing event with " + numHits + " hits...");

            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int sector = atofAdcBank.getByte("sector", hitIndex);
                int layer = atofAdcBank.getByte("layer", hitIndex);
                int component = atofAdcBank.getShort("component", hitIndex);
                int order = atofAdcBank.getByte("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                float time = atofAdcBank.getFloat("time", hitIndex);
                int pedestal = atofAdcBank.getShort("ped", hitIndex);
                double phi = 2 * Math.PI * component / 60;

                Hit hit = new Hit(sector, layer, component, order, adc, time, phi, pedestal);

                System.out.printf("Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                        sector, layer, component, order, adc, time, pedestal, phi);

                if (layer == 0) {
                    barHits.add(hit);
                } else if (layer >= 10 && layer <= 19) {
                    wedgeHits.add(hit);
                }
            }

            // Only proceed if there are at least 2 bar hits and 1 wedge hit
            if (barHits.size() >= 2 && wedgeHits.size() >= 1) {
                Hit barLeftHit = barHits.get(0);
                Hit barRightHit = barHits.get(1);

                double deltaT = barLeftHit.time - barRightHit.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeftHit.time - (zBar - BAR_LENGTH / 2) / VEFF,
                                       barRightHit.time - (zBar + BAR_LENGTH / 2) / VEFF);
                double phiBar = barLeftHit.phi;

                System.out.printf("Bar Z and TBar Calculation -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n",
                        zBar, tBar, phiBar);

                int clustersInEvent = 0;

                for (Hit wedgeHit : wedgeHits) {
                    double zWedge = (wedgeHit.component % N_WEDGE - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;
                    double deltaZ = Math.abs(zBar - zWedge);
                    double deltaPhi = Math.abs(phiBar - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    System.out.printf("Wedge Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                            wedgeHit.sector, wedgeHit.layer, wedgeHit.component, wedgeHit.order, wedgeHit.adc, wedgeHit.time, wedgeHit.pedestal, wedgeHit.phi);
                    System.out.printf("  Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", deltaZ, deltaPhi, deltaTime);

                    // Calculate proton angle if momentum components available (assuming hypothetical values here)
                    double px = 1.0, py = 1.0, pz = 1.0; // Replace with actual momentum data if available
                    double angle = Math.atan2(Math.sqrt(px * px + py * py), pz);
                    angleList.add(angle);

                    // Check if the hit combination meets threshold criteria for clustering
                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        clustersInEvent++;
                    }
                }

                // Only add clusters with a minimum size of 3 (2 bar hits + at least 1 wedge hit)
                if (clustersInEvent >= 3) {
                    clusterSizes.add((double) clustersInEvent);
                    eventIndices.add((double) eventIndex); // Ensure event index is added as Double
                }
            }

            eventIndex++;
        }
    }

    private void createPlots() {
        createScatterPlot("Delta Z vs Proton Angle", "Proton Angle (rad)", "Delta Z (mm)", angleList, deltaZList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", eventIndices, clusterSizes);
    }

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Double> xData, List<Double> yData) {
        XYSeries series = new XYSeries(title);
        for (int i = 0; i < xData.size(); i++) {
            series.add(xData.get(i), yData.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    static class Hit {
        int sector;
        int layer;
        int component;
        int order;
        int adc;
        double time;
        double phi;
        int pedestal;

        Hit(int sector, int layer, int component, int order, int adc, double time, double phi, int pedestal) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.time = time;
            this.phi = phi;
            this.pedestal = pedestal;
        }
    }
}

*/










/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm, threshold for delta Z
    private static final double PHI_THRESHOLD = 0.1; // rad, threshold for delta Phi
    private static final double TIME_THRESHOLD = 1.0; // ns, threshold for delta Time

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static List<Integer> eventIndices = new ArrayList<>();
    private static List<Double> angleList = new ArrayList<>(); // List to store angles for plotting

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

       ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN  demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN();
        demo.createPlots();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();
        int eventIndex = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            System.out.println("\nProcessing event with " + numHits + " hits...");

            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int sector = atofAdcBank.getByte("sector", hitIndex);
                int layer = atofAdcBank.getByte("layer", hitIndex);
                int component = atofAdcBank.getShort("component", hitIndex);
                int order = atofAdcBank.getByte("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                float time = atofAdcBank.getFloat("time", hitIndex);
                int pedestal = atofAdcBank.getShort("ped", hitIndex);
                double phi = 2 * Math.PI * component / 60;

                Hit hit = new Hit(sector, layer, component, order, adc, time, phi, pedestal);

                System.out.printf("Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                        sector, layer, component, order, adc, time, pedestal, phi);

                if (layer == 0) {
                    barHits.add(hit);
                } else if (layer >= 10 && layer <= 19) {
                    wedgeHits.add(hit);
                }
            }

            // Only proceed if there are at least 2 bar hits and 1 wedge hit
            if (barHits.size() >= 2 && wedgeHits.size() >= 1) {
                Hit barLeftHit = barHits.get(0);
                Hit barRightHit = barHits.get(1);

                double deltaT = barLeftHit.time - barRightHit.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeftHit.time - (zBar - BAR_LENGTH / 2) / VEFF,
                                       barRightHit.time - (zBar + BAR_LENGTH / 2) / VEFF);
                double phiBar = barLeftHit.phi;

                System.out.printf("Bar Z and TBar Calculation -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n",
                        zBar, tBar, phiBar);

                int clustersInEvent = 0;

                for (Hit wedgeHit : wedgeHits) {
                    double zWedge = (wedgeHit.component % N_WEDGE - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;
                    double deltaZ = Math.abs(zBar - zWedge);
                    double deltaPhi = Math.abs(phiBar - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    System.out.printf("Wedge Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                            wedgeHit.sector, wedgeHit.layer, wedgeHit.component, wedgeHit.order, wedgeHit.adc, wedgeHit.time, wedgeHit.pedestal, wedgeHit.phi);
                    System.out.printf("  Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", deltaZ, deltaPhi, deltaTime);

                    // Calculate proton angle if momentum components available (assuming hypothetical values here)
                    double px = 1.0, py = 1.0, pz = 1.0; // Replace with actual momentum data if available
                    double angle = Math.atan2(Math.sqrt(px * px + py * py), pz);
                    angleList.add(angle);

                    // Check if the hit combination meets threshold criteria for clustering
                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        clustersInEvent++;
                    }
                }

                // Only add clusters with a minimum size of 3 (2 bar hits + at least 1 wedge hit)
                if (clustersInEvent >= 3) {
                    clusterSizes.add(clustersInEvent);
                    eventIndices.add(eventIndex);
                }
            }

            eventIndex++;
        }
    }

    private void createPlots() {
        createScatterPlot("Delta Z vs Proton Angle", "Proton Angle (rad)", "Delta Z (mm)", angleList, deltaZList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", eventIndices, clusterSizes);
    }

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Double> xData, List<Double> yData) {
        XYSeries series = new XYSeries(title);
        for (int i = 0; i < xData.size(); i++) {
            series.add(xData.get(i), yData.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    static class Hit {
        int sector;
        int layer;
        int component;
        int order;
        int adc;
        double time;
        double phi;
        int pedestal;

        Hit(int sector, int layer, int component, int order, int adc, double time, double phi, int pedestal) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.time = time;
            this.phi = phi;
            this.pedestal = pedestal;
        }
    }
}

*/












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
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 200.0; // mm, threshold for delta Z
    private static final double PHI_THRESHOLD = 0.2; // rad, threshold for delta Phi
    private static final double TIME_THRESHOLD = 2.0; // ns, threshold for delta Time

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static List<Integer> eventIndices = new ArrayList<>();

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN  demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN();
        demo.createPlots();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();
        int eventIndex = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            System.out.println("\nProcessing event with " + numHits + " hits...");

            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int sector = atofAdcBank.getByte("sector", hitIndex);
                int layer = atofAdcBank.getByte("layer", hitIndex);
                int component = atofAdcBank.getShort("component", hitIndex);
                int order = atofAdcBank.getByte("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                float time = atofAdcBank.getFloat("time", hitIndex);
                int pedestal = atofAdcBank.getShort("ped", hitIndex);
                double phi = 2 * Math.PI * component / 60;

                Hit hit = new Hit(sector, layer, component, order, adc, time, phi, pedestal);

                System.out.printf("Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                        sector, layer, component, order, adc, time, pedestal, phi);

                if (layer == 0) {
                    barHits.add(hit);
                } else if (layer >= 10 && layer <= 19) {
                    wedgeHits.add(hit);
                }
            }

            // Only proceed if there are at least 2 bar hits and 1 wedge hit
            if (barHits.size() >= 2 && wedgeHits.size() >= 1) {
                Hit barLeftHit = barHits.get(0);
                Hit barRightHit = barHits.get(1);

                double deltaT = barLeftHit.time - barRightHit.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeftHit.time - (zBar - BAR_LENGTH / 2) / VEFF,
                                       barRightHit.time - (zBar + BAR_LENGTH / 2) / VEFF);
                double phiBar = barLeftHit.phi;

                System.out.printf("Bar Z and TBar Calculation -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n",
                        zBar, tBar, phiBar);

                int clustersInEvent = 0;

                for (Hit wedgeHit : wedgeHits) {
                    double zWedge = (wedgeHit.component % N_WEDGE - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;
                    double deltaZ = Math.abs(zBar - zWedge);
                    double deltaPhi = Math.abs(phiBar - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    System.out.printf("Wedge Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                            wedgeHit.sector, wedgeHit.layer, wedgeHit.component, wedgeHit.order, wedgeHit.adc, wedgeHit.time, wedgeHit.pedestal, wedgeHit.phi);
                    System.out.printf("  Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", deltaZ, deltaPhi, deltaTime);

                    // Check if the hit combination meets threshold criteria for clustering
                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        clustersInEvent++;
                    }
                }

                // Only add clusters with a minimum size of 3 (2 bar hits + at least 1 wedge hit)
                if (clustersInEvent >= 3) {
                    clusterSizes.add(clustersInEvent);
                    eventIndices.add(eventIndex);
                }
            }

            eventIndex++;
        }
    }

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes, eventIndices);
    }

    private void createHistogramPlot(String title, String xAxisLabel, List<Double> data) {
        HistogramDataset dataset = new HistogramDataset();
        double[] values = data.stream().mapToDouble(Double::doubleValue).toArray();
        dataset.addSeries(title, values, 50);

        JFreeChart histogram = ChartFactory.createHistogram(title, xAxisLabel, "Frequency", dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(histogram);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> data, List<Integer> indices) {
        XYSeries series = new XYSeries("Cluster Size");
        for (int i = 0; i < data.size(); i++) {
            series.add(indices.get(i), data.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    static class Hit {
        int sector;
        int layer;
        int component;
        int order;
        int adc;
        double time;
        double phi;
        int pedestal;

        Hit(int sector, int layer, int component, int order, int adc, double time, double phi, int pedestal) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.time = time;
            this.phi = phi;
            this.pedestal = pedestal;
        }
    }
}

*/







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
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm, threshold for delta Z
    private static final double PHI_THRESHOLD = 0.1; // rad, threshold for delta Phi
    private static final double TIME_THRESHOLD = 1.0; // ns, threshold for delta Time

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static List<Integer> eventIndices = new ArrayList<>();

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

       ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN();
        demo.createPlots();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();
        int eventIndex = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            System.out.println("\nProcessing event with " + numHits + " hits...");

            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int sector = atofAdcBank.getByte("sector", hitIndex);
                int layer = atofAdcBank.getByte("layer", hitIndex);
                int component = atofAdcBank.getShort("component", hitIndex);
                int order = atofAdcBank.getByte("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                float time = atofAdcBank.getFloat("time", hitIndex);
                int pedestal = atofAdcBank.getShort("ped", hitIndex);
                double phi = 2 * Math.PI * component / 60;

                Hit hit = new Hit(sector, layer, component, order, adc, time, phi, pedestal);

                System.out.printf("Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                        sector, layer, component, order, adc, time, pedestal, phi);

                if (layer == 0) {
                    barHits.add(hit);
                } else if (layer >= 10 && layer <= 19) {
                    wedgeHits.add(hit);
                }
            }

            // Process clusters if there are at least 2 bar hits and 1 wedge hit
            if (barHits.size() >= 2 && wedgeHits.size() >= 1) {
                Hit barLeftHit = barHits.get(0);
                Hit barRightHit = barHits.get(1);

                double deltaT = barLeftHit.time - barRightHit.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeftHit.time - (zBar - BAR_LENGTH / 2) / VEFF,
                                       barRightHit.time - (zBar + BAR_LENGTH / 2) / VEFF);
                double phiBar = barLeftHit.phi;

                System.out.printf("Bar Z and TBar Calculation -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n",
                        zBar, tBar, phiBar);

                int clustersInEvent = 0;

                for (Hit wedgeHit : wedgeHits) {
                    double zWedge = (wedgeHit.component % N_WEDGE - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;
                    double deltaZ = Math.abs(zBar - zWedge);
                    double deltaPhi = Math.abs(phiBar - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    System.out.printf("Wedge Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                            wedgeHit.sector, wedgeHit.layer, wedgeHit.component, wedgeHit.order, wedgeHit.adc, wedgeHit.time, wedgeHit.pedestal, wedgeHit.phi);
                    System.out.printf("  Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", deltaZ, deltaPhi, deltaTime);

                    // Check if this hit meets threshold criteria for clustering
                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        clustersInEvent++;
                    }
                }

                if (clustersInEvent > 0) {
                    clusterSizes.add(clustersInEvent);
                    eventIndices.add(eventIndex);
                }
            }

            eventIndex++;
        }
    }

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", clusterSizes, eventIndices);
    }

    private void createHistogramPlot(String title, String xAxisLabel, List<Double> data) {
        HistogramDataset dataset = new HistogramDataset();
        double[] values = data.stream().mapToDouble(Double::doubleValue).toArray();
        dataset.addSeries(title, values, 50);

        JFreeChart histogram = ChartFactory.createHistogram(title, xAxisLabel, "Frequency", dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(histogram);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> data, List<Integer> indices) {
        XYSeries series = new XYSeries("Cluster Size");
        for (int i = 0; i < data.size(); i++) {
            series.add(indices.get(i), data.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    static class Hit {
        int sector;
        int layer;
        int component;
        int order;
        int adc;
        double time;
        double phi;
        int pedestal;

        Hit(int sector, int layer, int component, int order, int adc, double time, double phi, int pedestal) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.time = time;
            this.phi = phi;
            this.pedestal = pedestal;
        }
    }
}

*/





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
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm, threshold for clustering
    private static final double PHI_THRESHOLD = 0.1; // rad, threshold for clustering
    private static final double TIME_THRESHOLD = 1.0; // ns, threshold for clustering

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static List<Integer> eventIndices = new ArrayList<>();

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN();
        demo.createPlots();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();
        int eventIndex = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            System.out.println("\nProcessing event with " + numHits + " hits...");

            Hit barLeftHit = null;
            Hit barRightHit = null;
            List<Hit> wedgeHits = new ArrayList<>();

            // Step 1: Identify bar and wedge hits within the event
            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int sector = atofAdcBank.getByte("sector", hitIndex);
                int layer = atofAdcBank.getByte("layer", hitIndex);
                int component = atofAdcBank.getShort("component", hitIndex);
                int order = atofAdcBank.getByte("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                float time = atofAdcBank.getFloat("time", hitIndex);
                int pedestal = atofAdcBank.getShort("ped", hitIndex);
                double phi = 2 * Math.PI * component / 60;

                System.out.printf("Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                        sector, layer, component, order, adc, time, pedestal, phi);

                if (layer == 0) { // Process bar hits
                    if (order == 0) {
                        barLeftHit = new Hit(sector, layer, component, order, adc, time, pedestal, phi);
                    } else if (order == 1) {
                        barRightHit = new Hit(sector, layer, component, order, adc, time, pedestal, phi);
                    }
                } else if (layer >= 10 && layer <= 19) { // Wedge hits
                    wedgeHits.add(new Hit(sector, layer, component, order, adc, time, pedestal, phi));
                }
            }

            int validClusterCount = 0;

            // Step 2: Calculate ZBar and TBar if both bar hits are present
            if (barLeftHit != null && barRightHit != null && barLeftHit.component == barRightHit.component) {
                double deltaT = barLeftHit.time - barRightHit.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeftHit.time - (zBar - BAR_LENGTH / 2) / VEFF,
                        barRightHit.time - (zBar + BAR_LENGTH / 2) / VEFF);

                System.out.printf("Bar Z and TBar Calculation -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n",
                        zBar, tBar, barLeftHit.phi);

                // Step 3: Match each wedge hit with the bar hits to form clusters
                for (Hit wedgeHit : wedgeHits) {
                    int wedgeIndexWithinBar = wedgeHit.component % N_WEDGE;
                    double zWedge = (wedgeIndexWithinBar - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;

                    double deltaZ = Math.abs(zBar - zWedge);
                    double deltaPhi = Math.abs(barLeftHit.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    System.out.printf("Wedge Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                            wedgeHit.sector, wedgeHit.layer, wedgeHit.component, wedgeHit.order, wedgeHit.adc, wedgeHit.time, wedgeHit.pedestal, wedgeHit.phi);
                    System.out.printf("  Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", deltaZ, deltaPhi, deltaTime);

                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    // Validate the cluster based on thresholds
                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        validClusterCount++;
                        int clusterSize = 3; // 2 bar hits + 1 wedge hit
                        clusterSizes.add(clusterSize);
                        eventIndices.add(eventIndex);
                        System.out.printf("Valid Cluster Formed -> Cluster Size: %d (2 bar hits + 1 wedge hit)\n", clusterSize);
                    }
                }
            } else {
                System.out.println("Invalid Cluster: Missing bar hits or wedge hits.");
            }

            // Print if no valid clusters were formed
            if (validClusterCount == 0) {
                System.out.println("No valid clusters formed in this event.");
            }

            eventIndex++;
        }
    }

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", eventIndices, clusterSizes);
    }

    private void createHistogramPlot(String title, String xAxisLabel, List<Double> data) {
        HistogramDataset dataset = new HistogramDataset();
        double[] values = data.stream().mapToDouble(Double::doubleValue).toArray();
        dataset.addSeries(title, values, 50);

        JFreeChart histogram = ChartFactory.createHistogram(title, xAxisLabel, "Frequency", dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(histogram);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData, List<Integer> yData) {
        XYSeries series = new XYSeries("Cluster Size");
        for (int i = 0; i < xData.size(); i++) {
            series.add(xData.get(i), yData.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    static class Hit {
        int sector;
        int layer;
        int component;
        int order;
        int adc;
        double time;
        int pedestal;
        double phi;

        Hit(int sector, int layer, int component, int order, int adc, double time, int pedestal, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.time = time;
            this.pedestal = pedestal;
            this.phi = phi;
        }
    }
}
*/






/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.ui.ApplicationFrame;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN extends ApplicationFrame {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm, threshold for clustering
    private static final double PHI_THRESHOLD = 0.1; // rad, threshold for clustering
    private static final double TIME_THRESHOLD = 1.0; // ns, threshold for clustering

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static List<Integer> eventIndices = new ArrayList<>();

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Cluster Analysis Plots");
        demo.createPlots();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();
        int eventIndex = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            System.out.println("\nProcessing event with " + numHits + " hits...");

            Hit barLeftHit = null;
            Hit barRightHit = null;
            List<Hit> wedgeHits = new ArrayList<>();

            // Step 1: Identify bar and wedge hits within the event
            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int sector = atofAdcBank.getByte("sector", hitIndex);
                int layer = atofAdcBank.getByte("layer", hitIndex);
                int component = atofAdcBank.getShort("component", hitIndex);
                int order = atofAdcBank.getByte("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                float time = atofAdcBank.getFloat("time", hitIndex);
                int pedestal = atofAdcBank.getShort("ped", hitIndex);
                double phi = 2 * Math.PI * component / 60;

                System.out.printf("Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                        sector, layer, component, order, adc, time, pedestal, phi);

                if (layer == 0) { // Process bar hits
                    if (order == 0) {
                        barLeftHit = new Hit(sector, layer, component, order, adc, time, pedestal, phi);
                    } else if (order == 1) {
                        barRightHit = new Hit(sector, layer, component, order, adc, time, pedestal, phi);
                    }
                } else if (layer >= 10 && layer <= 19) { // Wedge hits
                    wedgeHits.add(new Hit(sector, layer, component, order, adc, time, pedestal, phi));
                }
            }

            int validClusterCount = 0;

            // Step 2: Calculate ZBar and TBar if both bar hits are present
            if (barLeftHit != null && barRightHit != null && barLeftHit.component == barRightHit.component) {
                double deltaT = barLeftHit.time - barRightHit.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeftHit.time - (zBar - BAR_LENGTH / 2) / VEFF,
                        barRightHit.time - (zBar + BAR_LENGTH / 2) / VEFF);

                System.out.printf("Bar Z and TBar Calculation -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n",
                        zBar, tBar, barLeftHit.phi);

                // Step 3: Match each wedge hit with the bar hits to form clusters
                for (Hit wedgeHit : wedgeHits) {
                    int wedgeIndexWithinBar = wedgeHit.component % N_WEDGE;
                    double zWedge = (wedgeIndexWithinBar - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;

                    double deltaZ = Math.abs(zBar - zWedge);
                    double deltaPhi = Math.abs(barLeftHit.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    System.out.printf("Wedge Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                            wedgeHit.sector, wedgeHit.layer, wedgeHit.component, wedgeHit.order, wedgeHit.adc, wedgeHit.time, wedgeHit.pedestal, wedgeHit.phi);
                    System.out.printf("  Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", deltaZ, deltaPhi, deltaTime);

                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    // Validate the cluster based on thresholds
                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        validClusterCount++;
                        clusterSizes.add(3); // 2 bar hits + 1 wedge hit
                        eventIndices.add(eventIndex);
                        System.out.printf("Valid Cluster Formed -> Cluster Size: %d (2 bar hits + 1 wedge hit)\n", 3);
                    }
                }
            } else {
                System.out.println("Invalid Cluster: Missing bar hits or wedge hits.");
            }

            // Print no valid clusters if none were formed for the event
            if (validClusterCount == 0) {
                System.out.println("No valid clusters formed in this event.");
            }

            eventIndex++;
        }
    }

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", eventIndices, clusterSizes);
    }

    private void createHistogramPlot(String title, String xAxisLabel, List<Double> data) {
        HistogramDataset dataset = new HistogramDataset();
        double[] values = data.stream().mapToDouble(Double::doubleValue).toArray();
        dataset.addSeries(title, values, 50);

        JFreeChart histogram = ChartFactory.createHistogram(title, xAxisLabel, "Frequency", dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(histogram);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData, List<Integer> yData) {
        XYSeries series = new XYSeries("Cluster Size");
        for (int i = 0; i < xData.size(); i++) {
            series.add(xData.get(i), yData.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    static class Hit {
        int sector;
        int layer;
        int component;
        int order;
        int adc;
        double time;
        int pedestal;
        double phi;

        Hit(int sector, int layer, int component, int order, int adc, double time, int pedestal, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.time = time;
            this.pedestal = pedestal;
            this.phi = phi;
        }
    }
}


*/







/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.statistics.HistogramDataset;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN extends JFrame {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm, threshold for clustering
    private static final double PHI_THRESHOLD = 0.1; // rad, threshold for clustering
    private static final double TIME_THRESHOLD = 1.0; // ns, threshold for clustering

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static List<Integer> eventIndices = new ArrayList<>();

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Cluster Analysis Plots");
        demo.createPlots();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();
        int eventIndex = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            System.out.println("\nProcessing event with " + numHits + " hits...");

            Hit barLeftHit = null;
            Hit barRightHit = null;
            List<Hit> wedgeHits = new ArrayList<>();

            // Step 1: Identify bar and wedge hits within the event
            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int sector = atofAdcBank.getByte("sector", hitIndex);
                int layer = atofAdcBank.getByte("layer", hitIndex);
                int component = atofAdcBank.getShort("component", hitIndex);
                int order = atofAdcBank.getByte("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                float time = atofAdcBank.getFloat("time", hitIndex);
                int pedestal = atofAdcBank.getShort("ped", hitIndex);
                double phi = 2 * Math.PI * component / 60;

                System.out.printf("Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                        sector, layer, component, order, adc, time, pedestal, phi);

                if (layer == 0) { // Process bar hits
                    if (order == 0) {
                        barLeftHit = new Hit(sector, layer, component, order, adc, time, pedestal, phi);
                    } else if (order == 1) {
                        barRightHit = new Hit(sector, layer, component, order, adc, time, pedestal, phi);
                    }
                } else if (layer >= 10 && layer <= 19) { // Wedge hits
                    wedgeHits.add(new Hit(sector, layer, component, order, adc, time, pedestal, phi));
                }
            }

            int validClusterCount = 0;

            // Step 2: Calculate ZBar and TBar if both bar hits are present
            if (barLeftHit != null && barRightHit != null && barLeftHit.component == barRightHit.component) {
                double deltaT = barLeftHit.time - barRightHit.time;
                double zBar = VEFF * deltaT / 2;
                double tBar = Math.min(barLeftHit.time - (zBar - BAR_LENGTH / 2) / VEFF,
                        barRightHit.time - (zBar + BAR_LENGTH / 2) / VEFF);

                System.out.printf("Bar Z and TBar Calculation -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n",
                        zBar, tBar, barLeftHit.phi);

                // Step 3: Match each wedge hit with the bar hits to form clusters
                for (Hit wedgeHit : wedgeHits) {
                    int wedgeIndexWithinBar = wedgeHit.component % N_WEDGE;
                    double zWedge = (wedgeIndexWithinBar - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;

                    double deltaZ = Math.abs(zBar - zWedge);
                    double deltaPhi = Math.abs(barLeftHit.phi - wedgeHit.phi);
                    double deltaTime = Math.abs(tBar - wedgeHit.time);

                    System.out.printf("Wedge Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                            wedgeHit.sector, wedgeHit.layer, wedgeHit.component, wedgeHit.order, wedgeHit.adc, wedgeHit.time, wedgeHit.pedestal, wedgeHit.phi);
                    System.out.printf("  Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", deltaZ, deltaPhi, deltaTime);

                    deltaZList.add(deltaZ);
                    deltaPhiList.add(deltaPhi);
                    deltaTimeList.add(deltaTime);

                    // Validate the cluster based on thresholds
                    if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                        validClusterCount++;
                        clusterSizes.add(3); // 2 bar hits + 1 wedge hit
                        eventIndices.add(eventIndex);
                        System.out.printf("Valid Cluster Formed -> Cluster Size: %d (2 bar hits + 1 wedge hit)\n", 3);
                    }
                }
            } else {
                System.out.println("Invalid Cluster: Missing bar hits or wedge hits.");
            }

            // Print no valid clusters if none were formed for the event
            if (validClusterCount == 0) {
                System.out.println("No valid clusters formed in this event.");
            }

            eventIndex++;
        }
    }

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", eventIndices, clusterSizes);
    }

    private void createHistogramPlot(String title, String xAxisLabel, List<Double> data) {
        HistogramDataset dataset = new HistogramDataset();
        double[] values = data.stream().mapToDouble(Double::doubleValue).toArray();
        dataset.addSeries(title, values, 50);

        JFreeChart histogram = ChartFactory.createHistogram(title, xAxisLabel, "Frequency", dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(histogram);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData, List<Integer> yData) {
        XYSeries series = new XYSeries("Cluster Size");
        for (int i = 0; i < xData.size(); i++) {
            series.add(xData.get(i), yData.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    static class Hit {
        int sector;
        int layer;
        int component;
        int order;
        int adc;
        double time;
        int pedestal;
        double phi;

        Hit(int sector, int layer, int component, int order, int adc, double time, int pedestal, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.time = time;
            this.pedestal = pedestal;
            this.phi = phi;
        }
    }
}
*/














/* muc better but prints still some wrong info 
   /* 
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.statistics.HistogramDataset;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN extends JFrame {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm, threshold for clustering
    private static final double PHI_THRESHOLD = 0.1; // rad, threshold for clustering
    private static final double TIME_THRESHOLD = 1.0; // ns, threshold for clustering

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static List<Integer> eventIndices = new ArrayList<>();

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Cluster Analysis Plots");
        demo.createPlots();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();
        int eventIndex = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            int clustersInEvent = 0;
            boolean hasValidCluster = false;

            System.out.println("\nProcessing event with " + numHits + " hits...");

            Hit barLeftHit = null;
            Hit barRightHit = null;

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int sector = atofAdcBank.getByte("sector", hitIndex);
                int layer = atofAdcBank.getByte("layer", hitIndex);
                int component = atofAdcBank.getShort("component", hitIndex);
                int order = atofAdcBank.getByte("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                float time = atofAdcBank.getFloat("time", hitIndex);
                int pedestal = atofAdcBank.getShort("ped", hitIndex);
                double phi = 2 * Math.PI * component / 60;

                System.out.printf("Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                        sector, layer, component, order, adc, time, pedestal, phi);

                if (layer == 0) { // Process bar hits
                    if (order == 0) {
                        barLeftHit = new Hit(sector, layer, component, order, adc, time, phi);
                    } else {
                        barRightHit = new Hit(sector, layer, component, order, adc, time, phi);
                    }

                    if (barLeftHit != null && barRightHit != null && barLeftHit.component == barRightHit.component) {
                        double deltaT = barLeftHit.time - barRightHit.time;
                        double zBar = VEFF * deltaT / 2;
                        double tBar = Math.min(barLeftHit.time - (zBar - BAR_LENGTH / 2) / VEFF,
                                               barRightHit.time - (zBar + BAR_LENGTH / 2) / VEFF);

                        System.out.printf("Bar Z and TBar Calculation -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n",
                                zBar, tBar, phi);

                        barLeftHit = null;
                        barRightHit = null;

                        for (int wedgeIndex = 0; wedgeIndex < numHits; wedgeIndex++) {
                            int wedgeLayer = atofAdcBank.getByte("layer", wedgeIndex);
                            if (wedgeLayer >= 10 && wedgeLayer <= 19) {
                                int wedgeComponent = atofAdcBank.getShort("component", wedgeIndex);
                                double wedgeTime = atofAdcBank.getFloat("time", wedgeIndex);
                                double wedgePhi = phi;

                                int wedgeIndexWithinBar = wedgeComponent % N_WEDGE;
                                double zWedge = (wedgeIndexWithinBar - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;

                                double deltaZ = Math.abs(zBar - zWedge);
                                double deltaPhi = Math.abs(phi - wedgePhi);
                                double deltaTime = Math.abs(tBar - wedgeTime);

                                System.out.printf("Wedge Hit -> Sector: %d, Component: %d, ZWedge: %.2f mm, Time: %.2f ns, Phi: %.2f rad\n",
                                        sector, wedgeComponent, zWedge, wedgeTime, wedgePhi);
                                System.out.printf("  Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", deltaZ, deltaPhi, deltaTime);

                                deltaZList.add(deltaZ);
                                deltaPhiList.add(deltaPhi);
                                deltaTimeList.add(deltaTime);

                                boolean withinThreshold = deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD;
                                if (withinThreshold) {
                                    clustersInEvent++;
                                    hasValidCluster = true;
                                }
                            }
                        }
                    }
                }
            }

            // Determine if cluster is valid or invalid
            if (hasValidCluster && clustersInEvent >= 3) {
                System.out.println("Valid Cluster with size: " + clustersInEvent);
                clusterSizes.add(clustersInEvent);
                eventIndices.add(eventIndex);
            } else {
                System.out.println("Invalid Cluster");
            }

            eventIndex++;
        }
    }

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", eventIndices, clusterSizes);
    }

    private void createHistogramPlot(String title, String xAxisLabel, List<Double> data) {
        HistogramDataset dataset = new HistogramDataset();
        double[] values = data.stream().mapToDouble(Double::doubleValue).toArray();
        dataset.addSeries(title, values, 50);

        JFreeChart histogram = ChartFactory.createHistogram(title, xAxisLabel, "Frequency", dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(histogram);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData, List<Integer> yData) {
        XYSeries series = new XYSeries("Cluster Size");
        for (int i = 0; i < xData.size(); i++) {
            series.add(xData.get(i), yData.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    static class Hit {
        int sector;
        int layer;
        int component;
        int order;
        int adc;
        double time;
        double phi;

        Hit(int sector, int layer, int component, int order, int adc, double time, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.time = time;
            this.phi = phi;
        }
    }
}
*/


//plots correctly and prints wrong
/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.statistics.HistogramDataset;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN extends JFrame {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm, threshold for clustering
    private static final double PHI_THRESHOLD = 0.1; // rad, threshold for clustering
    private static final double TIME_THRESHOLD = 1.0; // ns, threshold for clustering

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();
    private static List<Integer> eventIndices = new ArrayList<>();

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

        ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Cluster Analysis Plots");
        demo.createPlots();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();
        int eventIndex = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            int clustersInEvent = 0;

            System.out.println("\nProcessing event with " + numHits + " hits...");

            Hit barLeftHit = null;
            Hit barRightHit = null;

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int sector = atofAdcBank.getByte("sector", hitIndex);
                int layer = atofAdcBank.getByte("layer", hitIndex);
                int component = atofAdcBank.getShort("component", hitIndex);
                int order = atofAdcBank.getByte("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                float time = atofAdcBank.getFloat("time", hitIndex);
                int pedestal = atofAdcBank.getShort("ped", hitIndex);
                double phi = 2 * Math.PI * component / 60;

                System.out.printf("Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                        sector, layer, component, order, adc, time, pedestal, phi);

                if (layer == 0) {
                    if (order == 0) {
                        barLeftHit = new Hit(sector, layer, component, order, adc, time, phi);
                    } else {
                        barRightHit = new Hit(sector, layer, component, order, adc, time, phi);
                    }

                    if (barLeftHit != null && barRightHit != null && barLeftHit.component == barRightHit.component) {
                        double deltaT = barLeftHit.time - barRightHit.time;
                        double zBar = VEFF * deltaT / 2;
                        double tBar = Math.min(barLeftHit.time - (zBar - BAR_LENGTH / 2) / VEFF,
                                               barRightHit.time - (zBar + BAR_LENGTH / 2) / VEFF);

                        System.out.printf("Bar Z and TBar Calculation -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n",
                                zBar, tBar, phi);

                        barLeftHit = null;
                        barRightHit = null;

                        for (int wedgeIndex = 0; wedgeIndex < numHits; wedgeIndex++) {
                            int wedgeLayer = atofAdcBank.getByte("layer", wedgeIndex);
                            if (wedgeLayer >= 10 && wedgeLayer <= 19) {
                                int wedgeComponent = atofAdcBank.getShort("component", wedgeIndex);
                                double wedgeTime = atofAdcBank.getFloat("time", wedgeIndex);
                                double wedgePhi = phi;

                                int wedgeIndexWithinBar = wedgeComponent % N_WEDGE;
                                double zWedge = (wedgeIndexWithinBar - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;

                                double deltaZ = Math.abs(zBar - zWedge);
                                double deltaPhi = Math.abs(phi - wedgePhi);
                                double deltaTime = Math.abs(tBar - wedgeTime);

                                System.out.printf("Wedge Hit -> Sector: %d, Component: %d, ZWedge: %.2f mm, Time: %.2f ns, Phi: %.2f rad\n",
                                        sector, wedgeComponent, zWedge, wedgeTime, wedgePhi);
                                System.out.printf("  Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", deltaZ, deltaPhi, deltaTime);

                                deltaZList.add(deltaZ);
                                deltaPhiList.add(deltaPhi);
                                deltaTimeList.add(deltaTime);

                                boolean isValidCluster = deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD;
                                if (isValidCluster) {
                                    clustersInEvent++;
                                }
                            }
                        }
                    }
                }
            }

            if (clustersInEvent < 3) {
                System.out.println("Invalid Cluster (size < 3)");
            } else {
                System.out.println("Valid Cluster with size: " + clustersInEvent);
                clusterSizes.add(clustersInEvent);
                eventIndices.add(eventIndex);
            }

            eventIndex++;
        }
    }

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", eventIndices, clusterSizes);
    }

    private void createHistogramPlot(String title, String xAxisLabel, List<Double> data) {
        HistogramDataset dataset = new HistogramDataset();
        double[] values = data.stream().mapToDouble(Double::doubleValue).toArray();
        dataset.addSeries(title, values, 50);

        JFreeChart histogram = ChartFactory.createHistogram(title, xAxisLabel, "Frequency", dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(histogram);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> xData, List<Integer> yData) {
        XYSeries series = new XYSeries("Cluster Size");
        for (int i = 0; i < xData.size(); i++) {
            series.add(xData.get(i), yData.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    static class Hit {
        int sector;
        int layer;
        int component;
        int order;
        int adc;
        double time;
        double phi;

        Hit(int sector, int layer, int component, int order, int adc, double time, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.time = time;
            this.phi = phi;
        }
    }
}

*/





































/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.ui.ApplicationFrame;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.statistics.HistogramDataset;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN  extends ApplicationFrame {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm, threshold for clustering
    private static final double PHI_THRESHOLD = 0.1; // rad, threshold for clustering
    private static final double TIME_THRESHOLD = 1.0; // ns, threshold for clustering

    private static List<Double> deltaZList = new ArrayList<>();
    private static List<Double> deltaPhiList = new ArrayList<>();
    private static List<Double> deltaTimeList = new ArrayList<>();
    private static List<Integer> clusterSizes = new ArrayList<>();

    public ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN(String title) {
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();

       ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN  demo = new ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN("Cluster Analysis Plots");
        demo.createPlots();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            int clustersInEvent = 0;

            Hit barLeftHit = null;
            Hit barRightHit = null;

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int sector = atofAdcBank.getByte("sector", hitIndex);
                int layer = atofAdcBank.getByte("layer", hitIndex);
                int component = atofAdcBank.getShort("component", hitIndex);
                int order = atofAdcBank.getByte("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                float time = atofAdcBank.getFloat("time", hitIndex);
                double phi = 2 * Math.PI * component / 60;

                if (layer == 0) { // Process bar hits only
                    if (order == 0) { // Left side of bar
                        barLeftHit = new Hit(sector, layer, component, order, adc, time, phi);
                    } else { // Right side of bar
                        barRightHit = new Hit(sector, layer, component, order, adc, time, phi);
                    }

                    if (barLeftHit != null && barRightHit != null && barLeftHit.component == barRightHit.component) {
                        double deltaT = barLeftHit.time - barRightHit.time;
                        double zBar = VEFF * deltaT / 2;
                        double tBar = Math.min(barLeftHit.time - (zBar - BAR_LENGTH / 2) / VEFF,
                                               barRightHit.time - (zBar + BAR_LENGTH / 2) / VEFF);

                        barLeftHit = null;
                        barRightHit = null;

                        for (int wedgeIndex = 0; wedgeIndex < numHits; wedgeIndex++) {
                            int wedgeLayer = atofAdcBank.getByte("layer", wedgeIndex);
                            if (wedgeLayer >= 10 && wedgeLayer <= 19) {
                                int wedgeComponent = atofAdcBank.getShort("component", wedgeIndex);
                                double wedgeTime = atofAdcBank.getFloat("time", wedgeIndex);
                                double wedgePhi = phi;

                                int wedgeIndexWithinBar = wedgeComponent % N_WEDGE;
                                double zWedge = (wedgeIndexWithinBar - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;

                                double deltaZ = Math.abs(zBar - zWedge);
                                double deltaPhi = Math.abs(phi - wedgePhi);
                                double deltaTime = Math.abs(tBar - wedgeTime);

                                boolean isValidCluster = deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD;
                                if (isValidCluster) {
                                    clustersInEvent++;
                                    deltaZList.add(deltaZ);
                                    deltaPhiList.add(deltaPhi);
                                    deltaTimeList.add(deltaTime);
                                }
                            }
                        }
                    }
                }
            }
            clusterSizes.add(clustersInEvent);
        }
    }

    private void createPlots() {
        createHistogramPlot("Delta Z Distribution", "Delta Z (mm)", deltaZList);
        createHistogramPlot("Delta Time Distribution", "Delta Time (ns)", deltaTimeList);
        createHistogramPlot("Delta Phi Distribution", "Delta Phi (rad)", deltaPhiList);
        createScatterPlot("Cluster Size vs Cluster Index", "Cluster Index", "Cluster Size", clusterSizes);
    }

    private void createHistogramPlot(String title, String xAxisLabel, List<Double> data) {
        HistogramDataset dataset = new HistogramDataset();
        double[] values = data.stream().mapToDouble(Double::doubleValue).toArray();
        dataset.addSeries(title, values, 50);

        JFreeChart histogram = ChartFactory.createHistogram(title, xAxisLabel, "Frequency", dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(histogram);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private void createScatterPlot(String title, String xAxisLabel, String yAxisLabel, List<Integer> data) {
        XYSeries series = new XYSeries("Cluster Size");
        for (int i = 0; i < data.size(); i++) {
            series.add(i, data.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart scatterPlot = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset,
                PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(scatterPlot);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    static class Hit {
        int sector;
        int layer;
        int component;
        int order;
        int adc;
        double time;
        double phi;

        Hit(int sector, int layer, int component, int order, int adc, double time, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.time = time;
            this.phi = phi;
        }
    }
}

*/








/*

package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 200.0; // mm, threshold for clustering
    private static final double PHI_THRESHOLD = 0.2; // rad, threshold for clustering
    private static final double TIME_THRESHOLD = 2.5; // ns, threshold for clustering

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();
        int clusterCount = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            System.out.println("\nProcessing event with " + numHits + " hits...");

            Hit barLeftHit = null;
            Hit barRightHit = null;

            // First pass: Identify and process bar hits (layer 0)
            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int sector = atofAdcBank.getByte("sector", hitIndex);
                int layer = atofAdcBank.getByte("layer", hitIndex);
                int component = atofAdcBank.getShort("component", hitIndex);
                int order = atofAdcBank.getByte("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                float time = atofAdcBank.getFloat("time", hitIndex);
                int ped = atofAdcBank.getShort("ped", hitIndex);
                double phi = 2 * Math.PI * component / 60;

                if (layer == 0) { // Process bar hits only
                    if (order == 0) { // Left side of bar
                        barLeftHit = new Hit(sector, layer, component, order, adc, time, ped, phi);
                        System.out.printf("Bar Left Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                                sector, layer, component, order, adc, time, ped, phi);
                    } else { // Right side of bar
                        barRightHit = new Hit(sector, layer, component, order, adc, time, ped, phi);
                        System.out.printf("Bar Right Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                                sector, layer, component, order, adc, time, ped, phi);
                    }

                    if (barLeftHit != null && barRightHit != null && barLeftHit.component == barRightHit.component) {
                        double deltaT = barLeftHit.time - barRightHit.time;
                        double zBar = VEFF * deltaT / 2;
                        double tBar = Math.min(barLeftHit.time - (zBar - BAR_LENGTH / 2) / VEFF, 
                                               barRightHit.time - (zBar + BAR_LENGTH / 2) / VEFF);

                        System.out.printf("Bar Z and TBar Calculation -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n", 
                                          zBar, tBar, phi);

                        barLeftHit = null;
                        barRightHit = null;

                        // Second pass: Identify wedge hits and calculate deltas for each wedge hit
                        for (int wedgeIndex = 0; wedgeIndex < numHits; wedgeIndex++) {
                            int wedgeSector = atofAdcBank.getByte("sector", wedgeIndex);
                            int wedgeLayer = atofAdcBank.getByte("layer", wedgeIndex);
                            int wedgeComponent = atofAdcBank.getShort("component", wedgeIndex);
                            int wedgeOrder = atofAdcBank.getByte("order", wedgeIndex);
                            int wedgeAdc = atofAdcBank.getInt("ADC", wedgeIndex);
                            float wedgeTime = atofAdcBank.getFloat("time", wedgeIndex);
                            int wedgePed = atofAdcBank.getShort("ped", wedgeIndex);
                            double wedgePhi = phi;

                            if (wedgeLayer >= 10 && wedgeLayer <= 19) { // Process wedge hits only
                                int wedgeIndexWithinBar = wedgeComponent % N_WEDGE;
                                double zWedge = (wedgeIndexWithinBar - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;

                                System.out.printf("Wedge Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, ZWedge: %.2f mm, Phi: %.2f rad\n",
                                                  wedgeSector, wedgeLayer, wedgeComponent, wedgeOrder, wedgeAdc, wedgeTime, wedgePed, zWedge, wedgePhi);

                                double deltaZ = Math.abs(zBar - zWedge);
                                double deltaPhi = Math.abs(phi - wedgePhi);
                                double deltaTime = Math.abs(tBar - wedgeTime);

                                System.out.printf("  Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", 
                                                  deltaZ, deltaPhi, deltaTime);

                                boolean isValidCluster = deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD;
                                if (isValidCluster) {
                                    clusterCount++;
                                    System.out.println("Valid Cluster Formed:");
                                    System.out.printf("  Bar Hits -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n", 
                                                      sector, layer, component, order, adc, time, ped, zBar, tBar, phi);
                                    System.out.printf("  Wedge Hit -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Z Position: %.2f mm, Phi: %.2f rad\n",
                                                      wedgeSector, wedgeLayer, wedgeComponent, wedgeOrder, wedgeAdc, wedgeTime, wedgePed, zWedge, wedgePhi);
                                } else {
                                    System.out.println("Invalid Cluster:");
                                    if (deltaZ >= Z_THRESHOLD) {
                                        System.out.printf("  Reason: Delta Z (%.2f mm) exceeds threshold of %.2f mm\n", deltaZ, Z_THRESHOLD);
                                    }
                                    if (deltaPhi >= PHI_THRESHOLD) {
                                        System.out.printf("  Reason: Delta Phi (%.2f rad) exceeds threshold of %.2f rad\n", deltaPhi, PHI_THRESHOLD);
                                    }
                                    if (deltaTime >= TIME_THRESHOLD) {
                                        System.out.printf("  Reason: Delta Time (%.2f ns) exceeds threshold of %.2f ns\n", deltaTime, TIME_THRESHOLD);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        System.out.printf("\nTotal Valid Clusters Formed Across All Events: %d\n", clusterCount);
    }

    static class Hit {
        int sector;
        int layer;
        int component;
        int order;
        int adc;
        double time;
        int ped;
        double phi;

        Hit(int sector, int layer, int component, int order, int adc, double time, int ped, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.time = time;
            this.ped = ped;
            this.phi = phi;
        }
    }
}





*/











// Sector , Layer, Component,ADC,  Time , Ped
/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm, threshold for clustering
    private static final double PHI_THRESHOLD = 0.1; // rad, threshold for clustering
    private static final double TIME_THRESHOLD = 1.0; // ns, threshold for clustering

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();
        int clusterCount = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            System.out.println("\nProcessing event with " + numHits + " hits...");

            Hit barLeftHit = null;
            Hit barRightHit = null;

            // First pass: Identify and process bar hits (layer 0)
            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int layer = atofAdcBank.getByte("layer", hitIndex);
                int sector = atofAdcBank.getByte("sector", hitIndex);
                int component = atofAdcBank.getShort("component", hitIndex);
                int order = atofAdcBank.getByte("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                float time = atofAdcBank.getFloat("time", hitIndex);
                double phi = 2 * Math.PI * component / 60;

                if (layer == 0) { // Process bar hits only
                    if (order == 0) { // Left side of bar
                        barLeftHit = new Hit(sector, layer, component, adc, time, phi);
                        System.out.printf("Bar Left Hit -> Sector: %d, Layer: %d, Component: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad\n",
                                sector, layer, component, adc, time, phi);
                    } else { // Right side of bar
                        barRightHit = new Hit(sector, layer, component, adc, time, phi);
                        System.out.printf("Bar Right Hit -> Sector: %d, Layer: %d, Component: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad\n",
                                sector, layer, component, adc, time, phi);
                    }

                    if (barLeftHit != null && barRightHit != null && barLeftHit.component == barRightHit.component) {
                        double deltaT = barLeftHit.time - barRightHit.time;
                        double zBar = VEFF * deltaT / 2;
                        double tBar = Math.min(barLeftHit.time - (zBar - BAR_LENGTH / 2) / VEFF, 
                                               barRightHit.time - (zBar + BAR_LENGTH / 2) / VEFF);

                        System.out.printf("Bar Z and TBar Calculation -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n", 
                                          zBar, tBar, phi);

                        barLeftHit = null;
                        barRightHit = null;

                        // Second pass: Identify wedge hits and calculate deltas for each wedge hit
                        for (int wedgeIndex = 0; wedgeIndex < numHits; wedgeIndex++) {
                            int wedgeLayer = atofAdcBank.getByte("layer", wedgeIndex);
                            if (wedgeLayer >= 10 && wedgeLayer <= 19) { // Process wedge hits only
                                int wedgeSector = atofAdcBank.getByte("sector", wedgeIndex);
                                int wedgeComponent = atofAdcBank.getShort("component", wedgeIndex); // Use component as shared identifier
                                double wedgeTime = atofAdcBank.getFloat("time", wedgeIndex);
                                double wedgePhi = phi; // Assuming the same phi as the bar
                                
                                // Calculate ZWedge based on wedge index within bar
                                int wedgeIndexWithinBar = wedgeComponent % N_WEDGE;
                                double zWedge = (wedgeIndexWithinBar - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;

                                System.out.printf("Wedge Hit -> Sector: %d, Layer: %d, Component: %d, ZWedge: %.2f mm, Time: %.2f ns, Phi: %.2f rad\n",
                                                  wedgeSector, wedgeLayer, wedgeComponent, zWedge, wedgeTime, wedgePhi);

                                double deltaZ = Math.abs(zBar - zWedge);
                                double deltaPhi = Math.abs(phi - wedgePhi);
                                double deltaTime = Math.abs(tBar - wedgeTime);

                                // Print delta values for each bar-wedge combination
                                System.out.printf("  Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", 
                                                  deltaZ, deltaPhi, deltaTime);

                                // Check if this bar-wedge combination forms a valid cluster
                                boolean isValidCluster = deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD;
                                if (isValidCluster) {
                                    clusterCount++;
                                    System.out.println("Valid Cluster Formed:");
                                    System.out.printf("  Bar Hits -> Sector: %d, Layer: %d, Component: %d, ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n", 
                                                      sector, layer, component, zBar, tBar, phi);
                                    System.out.printf("  Wedge Hit -> Sector: %d, Layer: %d, Component: %d, Z Position: %.2f mm, Time: %.2f ns, Phi: %.2f rad\n",
                                                      wedgeSector, wedgeLayer, wedgeComponent, zWedge, wedgeTime, wedgePhi);
                                } else {
                                    System.out.println("Invalid Cluster:");
                                    if (deltaZ >= Z_THRESHOLD) {
                                        System.out.printf("  Reason: Delta Z (%.2f mm) exceeds threshold of %.2f mm\n", deltaZ, Z_THRESHOLD);
                                    }
                                    if (deltaPhi >= PHI_THRESHOLD) {
                                        System.out.printf("  Reason: Delta Phi (%.2f rad) exceeds threshold of %.2f rad\n", deltaPhi, PHI_THRESHOLD);
                                    }
                                    if (deltaTime >= TIME_THRESHOLD) {
                                        System.out.printf("  Reason: Delta Time (%.2f ns) exceeds threshold of %.2f ns\n", deltaTime, TIME_THRESHOLD);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        System.out.printf("\nTotal Valid Clusters Formed Across All Events: %d\n", clusterCount);
    }

    static class Hit {
        int sector;
        int layer;
        int component;
        int adc;
        double time;
        double phi;

        Hit(int sector, int layer, int component, int adc, double time, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.adc = adc;
            this.time = time;
            this.phi = phi;
        }
    }
}
*/


















//Sector and  component 

/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm, threshold for clustering
    private static final double PHI_THRESHOLD = 0.1; // rad, threshold for clustering
    private static final double TIME_THRESHOLD = 1.0; // ns, threshold for clustering

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();
        int clusterCount = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            System.out.println("\nProcessing event with " + numHits + " hits...");

            Hit barLeftHit = null;
            Hit barRightHit = null;

            // First pass: Identify and process bar hits (layer 0)
            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int layer = atofAdcBank.getByte("layer", hitIndex);
                int sector = atofAdcBank.getByte("sector", hitIndex);
                int component = atofAdcBank.getShort("component", hitIndex);
                int order = atofAdcBank.getByte("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                float time = atofAdcBank.getFloat("time", hitIndex);
                double phi = 2 * Math.PI * component / 60;

                if (layer == 0) { // Process bar hits only
                    if (order == 0) { // Left side of bar
                        barLeftHit = new Hit(sector, layer, component, adc, time, phi);
                        System.out.printf("Bar Left Hit -> Sector: %d, Component: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad\n",
                                sector, component, adc, time, phi);
                    } else { // Right side of bar
                        barRightHit = new Hit(sector, layer, component, adc, time, phi);
                        System.out.printf("Bar Right Hit -> Sector: %d, Component: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad\n",
                                sector, component, adc, time, phi);
                    }

                    if (barLeftHit != null && barRightHit != null && barLeftHit.component == barRightHit.component) {
                        double deltaT = barLeftHit.time - barRightHit.time;
                        double zBar = VEFF * deltaT / 2;
                        double tBar = Math.min(barLeftHit.time - (zBar - BAR_LENGTH / 2) / VEFF, 
                                               barRightHit.time - (zBar + BAR_LENGTH / 2) / VEFF);

                        System.out.printf("Bar Z and TBar Calculation -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n", 
                                          zBar, tBar, phi);

                        barLeftHit = null;
                        barRightHit = null;

                        // Second pass: Identify wedge hits and calculate deltas for each wedge hit
                        for (int wedgeIndex = 0; wedgeIndex < numHits; wedgeIndex++) {
                            int wedgeLayer = atofAdcBank.getByte("layer", wedgeIndex);
                            if (wedgeLayer >= 10 && wedgeLayer <= 19) { // Process wedge hits only
                                int wedgeSector = atofAdcBank.getByte("sector", wedgeIndex);
                                int wedgeComponent = atofAdcBank.getShort("component", wedgeIndex); // Use component as shared identifier
                                double wedgeTime = atofAdcBank.getFloat("time", wedgeIndex);
                                double wedgePhi = phi; // Assuming the same phi as the bar
                                
                                // Calculate ZWedge based on wedge index within bar
                                int wedgeIndexWithinBar = wedgeComponent % N_WEDGE;
                                double zWedge = (wedgeIndexWithinBar - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;

                                System.out.printf("Wedge Hit -> Sector: %d, Component: %d, Layer: %d, ZWedge: %.2f mm, Time: %.2f ns, Phi: %.2f rad\n",
                                                  wedgeSector, wedgeComponent, wedgeLayer, zWedge, wedgeTime, wedgePhi);

                                double deltaZ = Math.abs(zBar - zWedge);
                                double deltaPhi = Math.abs(phi - wedgePhi);
                                double deltaTime = Math.abs(tBar - wedgeTime);

                                // Print delta values for each bar-wedge combination
                                System.out.printf("  Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", 
                                                  deltaZ, deltaPhi, deltaTime);

                                // Check if this bar-wedge combination forms a valid cluster
                                boolean isValidCluster = deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD;
                                if (isValidCluster) {
                                    clusterCount++;
                                    System.out.println("Valid Cluster Formed:");
                                    System.out.printf("  Bar Hits -> Sector: %d, Component: %d, ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n", 
                                                      sector, component, zBar, tBar, phi);
                                    System.out.printf("  Wedge Hit -> Sector: %d, Component: %d, Z Position: %.2f mm, Time: %.2f ns, Phi: %.2f rad\n",
                                                      wedgeSector, wedgeComponent, zWedge, wedgeTime, wedgePhi);
                                } else {
                                    System.out.println("Invalid Cluster:");
                                    if (deltaZ >= Z_THRESHOLD) {
                                        System.out.printf("  Reason: Delta Z (%.2f mm) exceeds threshold of %.2f mm\n", deltaZ, Z_THRESHOLD);
                                    }
                                    if (deltaPhi >= PHI_THRESHOLD) {
                                        System.out.printf("  Reason: Delta Phi (%.2f rad) exceeds threshold of %.2f rad\n", deltaPhi, PHI_THRESHOLD);
                                    }
                                    if (deltaTime >= TIME_THRESHOLD) {
                                        System.out.printf("  Reason: Delta Time (%.2f ns) exceeds threshold of %.2f ns\n", deltaTime, TIME_THRESHOLD);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        System.out.printf("\nTotal Valid Clusters Formed Across All Events: %d\n", clusterCount);
    }

    static class Hit {
        int sector;
        int layer;
        int component;
        int adc;
        double time;
        double phi;

        Hit(int sector, int layer, int component, int adc, double time, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.adc = adc;
            this.time = time;
            this.phi = phi;
        }
    }
}


*/













/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final int N_WEDGE = 10; // Total number of wedges per bar
    private static final double Z_THRESHOLD = 150.0; // mm, threshold for clustering
    private static final double PHI_THRESHOLD = 0.1; // rad, threshold for clustering
    private static final double TIME_THRESHOLD = 1.0; // ns, threshold for clustering

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();
        int clusterCount = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            System.out.println("\nProcessing event with " + numHits + " hits...");

            Hit barLeftHit = null;
            Hit barRightHit = null;

            // First pass: Identify and process bar hits (layer 0)
            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int layer = atofAdcBank.getByte("layer", hitIndex);
                if (layer == 0) { // Process bar hits only
                    int sector = atofAdcBank.getByte("sector", hitIndex);
                    int component = atofAdcBank.getShort("component", hitIndex);
                    int order = atofAdcBank.getByte("order", hitIndex);
                    int adc = atofAdcBank.getInt("ADC", hitIndex);
                    float time = atofAdcBank.getFloat("time", hitIndex);
                    double phiBar = 2 * Math.PI * component / 60;

                    if (order == 0) { // Left side
                        barLeftHit = new Hit(sector, layer, component, adc, time, phiBar);
                        System.out.printf("Bar Left Hit -> Sector: %d, Component: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad\n",
                                sector, component, adc, time, phiBar);
                    } else { // Right side
                        barRightHit = new Hit(sector, layer, component, adc, time, phiBar);
                        System.out.printf("Bar Right Hit -> Sector: %d, Component: %d, ADC: %d, Time: %.2f ns, Phi: %.2f rad\n",
                                sector, component, adc, time, phiBar);
                    }

                    if (barLeftHit != null && barRightHit != null && barLeftHit.component == barRightHit.component) {
                        double deltaT = barLeftHit.time - barRightHit.time;
                        double zBar = VEFF * deltaT / 2;
                        double tBar = Math.min(barLeftHit.time - (zBar - BAR_LENGTH / 2) / VEFF, 
                                               barRightHit.time - (zBar + BAR_LENGTH / 2) / VEFF);

                        System.out.printf("Bar Z and TBar Calculation -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n", 
                                          zBar, tBar, phiBar);

                        barLeftHit = null;
                        barRightHit = null;

                        // Second pass: Identify wedge hits and calculate deltas for each wedge hit
                        for (int wedgeIndex = 0; wedgeIndex < numHits; wedgeIndex++) {
                            int wedgeLayer = atofAdcBank.getByte("layer", wedgeIndex);
                            if (wedgeLayer >= 10 && wedgeLayer <= 19) { // Process wedge hits only
                                int wedgeSector = atofAdcBank.getByte("sector", wedgeIndex);
                                int wedgeIndexWithinBar = atofAdcBank.getShort("component", wedgeIndex) % N_WEDGE; // Ensure wedge index is 0-9
                                double wedgeTime = atofAdcBank.getFloat("time", wedgeIndex);
                                double wedgePhi = phiBar; // Assuming the same phi as the bar
                                
                                // Corrected Z_Wedge calculation using wedgeIndexWithinBar
                                double zWedge = (wedgeIndexWithinBar - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;

                                System.out.printf("Wedge Hit -> Sector: %d, Index: %d, Layer: %d, ZWedge: %.2f mm, Time: %.2f ns, Phi: %.2f rad\n",
                                                  wedgeSector, wedgeIndexWithinBar, wedgeLayer, zWedge, wedgeTime, wedgePhi);

                                double deltaZ = Math.abs(zBar - zWedge);
                                double deltaPhi = Math.abs(phiBar - wedgePhi);
                                double deltaTime = Math.abs(tBar - wedgeTime);

                                // Print delta values for each bar-wedge combination
                                System.out.printf("  Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", 
                                                  deltaZ, deltaPhi, deltaTime);

                                // Check if this bar-wedge combination forms a valid cluster
                                if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                                    clusterCount++;
                                    System.out.println("Valid Cluster Formed:");
                                    System.out.printf("  Bar Hits -> Sector: %d, Component: %d, ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n", 
                                                      sector, component, zBar, tBar, phiBar);
                                    System.out.printf("  Wedge Hit -> Sector: %d, Index: %d, Z Position: %.2f mm, Time: %.2f ns, Phi: %.2f rad\n",
                                                      wedgeSector, wedgeIndexWithinBar, zWedge, wedgeTime, wedgePhi);
                                }
                            }
                        }
                    }
                }
            }
        }
        System.out.printf("\nTotal Valid Clusters Formed Across All Events: %d\n", clusterCount);
    }

    static class Hit {
        int sector;
        int layer;
        int component;
        int adc;
        double time;
        double phi;

        Hit(int sector, int layer, int component, int adc, double time, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.adc = adc;
            this.time = time;
            this.phi = phi;
        }
    }
}
*/










/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final double Z_THRESHOLD = 150.0; // mm, threshold for clustering
    private static final double PHI_THRESHOLD = 0.1; // rad, threshold for clustering
    private static final double TIME_THRESHOLD = 1.0; // ns, threshold for clustering

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();
        int clusterCount = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            System.out.println("\nProcessing event with " + numHits + " hits...");

            Hit barLeftHit = null;
            Hit barRightHit = null;

            // First pass: Identify and process bar hits (layer 0)
            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int layer = atofAdcBank.getByte("layer", hitIndex);
                if (layer == 0) { // Process bar hits only
                    int sector = atofAdcBank.getByte("sector", hitIndex);
                    int component = atofAdcBank.getShort("component", hitIndex);
                    int order = atofAdcBank.getByte("order", hitIndex);
                    int adc = atofAdcBank.getInt("ADC", hitIndex);
                    float time = atofAdcBank.getFloat("time", hitIndex);
                    double phiBar = 2 * Math.PI * component / 60;

                    if (order == 0) { // Left side
                        barLeftHit = new Hit(sector, layer, component, adc, time, phiBar);
                        System.out.printf("Bar Left Hit -> Sector: %d, Component: %d, ADC: %d, Time: %.2f ns\n",
                                sector, component, adc, time);
                    } else { // Right side
                        barRightHit = new Hit(sector, layer, component, adc, time, phiBar);
                        System.out.printf("Bar Right Hit -> Sector: %d, Component: %d, ADC: %d, Time: %.2f ns\n",
                                sector, component, adc, time);
                    }

                    if (barLeftHit != null && barRightHit != null && barLeftHit.component == barRightHit.component) {
                        double deltaT = barLeftHit.time - barRightHit.time;
                        double zBar = VEFF * deltaT / 2;
                        double tBar = Math.min(barLeftHit.time - (zBar - BAR_LENGTH / 2) / VEFF, 
                                               barRightHit.time - (zBar + BAR_LENGTH / 2) / VEFF);

                        System.out.printf("Bar Z and TBar Calculation -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n", 
                                          zBar, tBar, phiBar);

                        barLeftHit = null;
                        barRightHit = null;

                        // Second pass: Identify wedge hits and check for clustering
                        for (int wedgeIndex = 0; wedgeIndex < numHits; wedgeIndex++) {
                            int wedgeLayer = atofAdcBank.getByte("layer", wedgeIndex);
                            if (wedgeLayer >= 10 && wedgeLayer <= 19) { // Process wedge hits only
                                int wedgeSector = atofAdcBank.getByte("sector", wedgeIndex);
                                int wedgeComponent = atofAdcBank.getShort("component", wedgeIndex);
                                double wedgeTime = atofAdcBank.getFloat("time", wedgeIndex);
                                double zWedge = (wedgeComponent - 5) * WEDGE_SPACING;
                                double wedgePhi = phiBar;

                                System.out.printf("Wedge Hit -> Sector: %d, Component: %d, Layer: %d, Time: %.2f ns\n",
                                                  wedgeSector, wedgeComponent, wedgeLayer, wedgeTime);

                                double deltaZ = Math.abs(zBar - zWedge);
                                double deltaPhi = Math.abs(phiBar - wedgePhi);
                                double deltaTime = Math.abs(tBar - wedgeTime);

                                if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                                    clusterCount++;
                                    System.out.println("Valid Cluster Formed:");
                                    System.out.printf("  Bar Hits -> Sector: %d, Component: %d, ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n", 
                                                      sector, component, zBar, tBar, phiBar);
                                    System.out.printf("  Wedge Hit -> Sector: %d, Component: %d, Z Position: %.2f mm, Time: %.2f ns, Phi: %.2f rad\n",
                                                      wedgeSector, wedgeComponent, zWedge, wedgeTime, wedgePhi);
                                    System.out.printf("  Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", 
                                                      deltaZ, deltaPhi, deltaTime);
                                }
                            }
                        }
                    }
                }
            }
        }
        System.out.printf("\nTotal Valid Clusters Formed Across All Events: %d\n", clusterCount);
    }

    static class Hit {
        int sector;
        int layer;
        int component;
        int adc;
        double time;
        double phi;

        Hit(int sector, int layer, int component, int adc, double time, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.adc = adc;
            this.time = time;
            this.phi = phi;
        }
    }
}


*/












/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final double Z_THRESHOLD = 150.0; // mm, threshold for clustering
    private static final double PHI_THRESHOLD = 0.1; // rad, threshold for clustering
    private static final double TIME_THRESHOLD = 1.0; // ns, threshold for clustering

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();
        int clusterCount = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            System.out.println("\nProcessing event with " + numHits + " hits...");

            Hit barLeftHit = null;
            Hit barRightHit = null;

            // First pass: Identify and process bar hits (layer 0)
            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int layer = atofAdcBank.getByte("layer", hitIndex);
                if (layer == 0) { // Process bar hits only
                    int sector = atofAdcBank.getByte("sector", hitIndex);
                    int component = atofAdcBank.getShort("component", hitIndex);
                    int order = atofAdcBank.getByte("order", hitIndex);
                    int adc = atofAdcBank.getInt("ADC", hitIndex);
                    float time = atofAdcBank.getFloat("time", hitIndex);
                    double phiBar = 2 * Math.PI * component / 60;

                    if (order == 0) { // Left side
                        barLeftHit = new Hit(sector, layer, component, adc, time, phiBar);
                        System.out.printf("Bar Left Hit -> Sector: %d, Component: %d, ADC: %d, Time: %.2f ns\n",
                                sector, component, adc, time);
                    } else { // Right side
                        barRightHit = new Hit(sector, layer, component, adc, time, phiBar);
                        System.out.printf("Bar Right Hit -> Sector: %d, Component: %d, ADC: %d, Time: %.2f ns\n",
                                sector, component, adc, time);
                    }

                    if (barLeftHit != null && barRightHit != null && barLeftHit.component == barRightHit.component) {
                        double deltaT = barLeftHit.time - barRightHit.time;
                        double zBar = VEFF * deltaT / 2;
                        double tBar = Math.min(barLeftHit.time - (zBar - BAR_LENGTH / 2) / VEFF, 
                                               barRightHit.time - (zBar + BAR_LENGTH / 2) / VEFF);

                        System.out.printf("Bar Z and TBar Calculation -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n", 
                                          zBar, tBar, phiBar);

                        barLeftHit = null;
                        barRightHit = null;

                        // Second pass: Identify wedge hits and check for clustering
                        for (int wedgeIndex = 0; wedgeIndex < numHits; wedgeIndex++) {
                            int wedgeLayer = atofAdcBank.getByte("layer", wedgeIndex);
                            if (wedgeLayer >= 10 && wedgeLayer <= 19) { // Process wedge hits only
                                int wedgeSector = atofAdcBank.getByte("sector", wedgeIndex);
                                int wedgeComponent = atofAdcBank.getShort("component", wedgeIndex);
                                double wedgeTime = atofAdcBank.getFloat("time", wedgeIndex);
                                double zWedge = (wedgeComponent - 5) * WEDGE_SPACING;
                                double wedgePhi = phiBar;

                                System.out.printf("Wedge Hit -> Sector: %d, Component: %d, Layer: %d, Time: %.2f ns\n",
                                                  wedgeSector, wedgeComponent, wedgeLayer, wedgeTime);

                                double deltaZ = Math.abs(zBar - zWedge);
                                double deltaPhi = Math.abs(phiBar - wedgePhi);
                                double deltaTime = Math.abs(tBar - wedgeTime);

                                if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                                    clusterCount++;
                                    System.out.println("Valid Cluster Formed:");
                                    System.out.printf("  Bar Hits -> Sector: %d, Component: %d, ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n", 
                                                      sector, component, zBar, tBar, phiBar);
                                    System.out.printf("  Wedge Hit -> Sector: %d, Component: %d, Z Position: %.2f mm, Time: %.2f ns, Phi: %.2f rad\n",
                                                      wedgeSector, wedgeComponent, zWedge, wedgeTime, wedgePhi);
                                    System.out.printf("  Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", 
                                                      deltaZ, deltaPhi, deltaTime);
                                }
                            }
                        }
                    }
                }
            }
        }
        System.out.printf("\nTotal Valid Clusters Formed Across All Events: %d\n", clusterCount);
    }

    static class Hit {
        int sector;
        int layer;
        int component;
        int adc;
        double time;
        double phi;

        Hit(int sector, int layer, int component, int adc, double time, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.adc = adc;
            this.time = time;
            this.phi = phi;
        }
    }
}

*/




















/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final double Z_THRESHOLD = 150.0; // mm, threshold for clustering
    private static final double PHI_THRESHOLD = 0.1; // rad, threshold for clustering
    private static final double TIME_THRESHOLD = 1.0; // ns, threshold for clustering

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();
        int clusterCount = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            System.out.println("\nProcessing event with " + numHits + " hits...");

            // Track bar hits
            Hit barLeftHit = null;
            Hit barRightHit = null;
            Hit wedgeHit = null;

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int sector = atofAdcBank.getByte("sector", hitIndex);
                int layer = atofAdcBank.getByte("layer", hitIndex);
                int component = atofAdcBank.getShort("component", hitIndex);
                int order = atofAdcBank.getByte("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                float time = atofAdcBank.getFloat("time", hitIndex);
                int pedestal = atofAdcBank.getShort("ped", hitIndex);

                if (layer == 0) { // Bar hit
                    double phiBar = 2 * Math.PI * component / 60;

                    if (order == 0) { // Left side
                        barLeftHit = new Hit(sector, layer, component, adc, time, phiBar);
                        System.out.printf("Bar Left Hit -> Sector: %d, Layer: %d, Component: %d, ADC: %d, Time: %.2f ns, Pedestal: %d\n",
                                sector, layer, component, adc, time, pedestal);
                    } else { // Right side
                        barRightHit = new Hit(sector, layer, component, adc, time, phiBar);
                        System.out.printf("Bar Right Hit -> Sector: %d, Layer: %d, Component: %d, ADC: %d, Time: %.2f ns, Pedestal: %d\n",
                                sector, layer, component, adc, time, pedestal);
                    }

                    if (barLeftHit != null && barRightHit != null && barLeftHit.component == barRightHit.component) {
                        // Calculate ZBar and TBar based on time difference
                        double deltaT = barLeftHit.time - barRightHit.time;
                        double zBar = VEFF * deltaT / 2;

                        // Calculate zLeft and zRight based on the new zBar
                        double zLeft = zBar - BAR_LENGTH / 2;
                        double zRight = zBar + BAR_LENGTH / 2;

                        // Calculate tBarLeft and tBarRight based on zLeft and zRight
                        double tBarLeft = barLeftHit.time - (zLeft / VEFF);
                        double tBarRight = barRightHit.time - (zRight / VEFF);
                        double tBar = Math.min(tBarLeft, tBarRight);

                        System.out.printf("Bar Z and TBar Calculation -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n", zBar, tBar, barLeftHit.phi);

                        // Reset barLeftHit and barRightHit for the next component
                        barLeftHit = null;
                        barRightHit = null;

                        // Now check for clustering with each wedge hit
                        for (int wedgeIndex = 0; wedgeIndex < numHits; wedgeIndex++) {
                            int wedgeLayer = atofAdcBank.getByte("layer", wedgeIndex);
                            if (wedgeLayer >= 10 && wedgeLayer <= 19) { // Wedge hit
                                int wedgeSector = atofAdcBank.getByte("sector", wedgeIndex);
                                int wedgeComponent = atofAdcBank.getShort("component", wedgeIndex) % 10;
                                double wedgeTime = atofAdcBank.getFloat("time", wedgeIndex);
                                double wedgePhi = phiBar; // Assume PhiBar applies to the wedge
                                double zWedge = (wedgeComponent - 5) * WEDGE_SPACING;

                                double deltaZ = Math.abs(zBar - zWedge);
                                double deltaPhi = Math.abs(phiBar - wedgePhi);
                                double deltaTime = Math.abs(tBar - wedgeTime);

                                if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                                    clusterCount++;
                                    System.out.println("Valid Cluster Formed:");
                                    System.out.printf("  Bar Hits -> Sector: %d, Component: %d\n", sector, component);
                                    System.out.printf("    Bar Hit -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n",
                                            zBar, tBar, phiBar);
                                    System.out.printf("  Wedge Hit -> Sector: %d, Component: %d, Z Position: %.2f mm, Time: %.2f ns, Phi: %.2f rad\n",
                                            wedgeSector, wedgeComponent, zWedge, wedgeTime, wedgePhi);
                                    System.out.printf("  Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", deltaZ, deltaPhi, deltaTime);
                                }
                            }
                        }
                    }
                }
            }
        }
        System.out.printf("\nTotal Valid Clusters Formed Across All Events: %d\n", clusterCount);
    }

    static class Hit {
        int sector;
        int layer;
        int component;
        int adc;
        double time;
        double phi;
        double zPosition;

        Hit(int sector, int layer, int component, int adc, double time, double phi) {
            this(sector, layer, component, adc, time, phi, 0);
        }

        Hit(int sector, int layer, int component, int adc, double time, double phi, double zPosition) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.adc = adc;
            this.time = time;
            this.phi = phi;
            this.zPosition = zPosition;
        }
    }
}


*/












/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import java.util.HashMap;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final double Z_THRESHOLD = 200.0; // mm, threshold for clustering
    private static final double PHI_THRESHOLD = 0.2; // rad, threshold for clustering
    private static final double TIME_THRESHOLD = 2.0; // ns, threshold for clustering

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();
        int clusterCount = 0;
        Map<Integer, Integer> clusterSizeCount = new HashMap<>();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            System.out.println("\nProcessing event with " + numHits + " hits...");

            boolean leftHitFound = false;
            boolean rightHitFound = false;
            double tBarLeft = 0.0;
            double tBarRight = 0.0;
            double zBar = 0.0;
            double tBar = 0.0;
            double phiBar = 0.0;
            int lastComponent = -1;

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int sector = atofAdcBank.getByte("sector", hitIndex);
                int layer = atofAdcBank.getByte("layer", hitIndex);
                int component = atofAdcBank.getShort("component", hitIndex);
                int order = atofAdcBank.getByte("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                float time = atofAdcBank.getFloat("time", hitIndex);
                int pedestal = atofAdcBank.getShort("ped", hitIndex);

                if (layer == 0) { // Bar hit
                    phiBar = 2 * Math.PI * component / 60; // Phi calculation for the bar and all wedges in this bar
                    if (order == 0) { // Left side
                        double zLeft = zBar - (BAR_LENGTH / 2); // Z position for left side
                        tBarLeft = time - (zLeft / VEFF); // Calculate TBar from Left
                        leftHitFound = true;
                        System.out.printf("Bar Hit -> Sector: %d, Layer: %d, Component: %d, ADC: %d, Time (Left): %.2f ns, Pedestal: %d\n",
                                          sector, layer, component, adc, time, pedestal);
                    } else { // Right side
                        double zRight = zBar + (BAR_LENGTH / 2); // Z position for right side
                        tBarRight = time - (zRight / VEFF); // Calculate TBar from Right
                        rightHitFound = true;
                        System.out.printf("Bar Hit -> Sector: %d, Layer: %d, Component: %d, ADC: %d, Time (Right): %.2f ns, Pedestal: %d\n",
                                          sector, layer, component, adc, time, pedestal);
                    }

                    if (leftHitFound && rightHitFound && lastComponent != component) { // Calculate ZBar and TBar after both hits are found
                        double deltaT = Math.abs(tBarLeft - tBarRight);
                        zBar = VEFF * deltaT / 2; // Calculate ZBar
                        tBar = Math.min(tBarLeft, tBarRight); // TBar using the minimum of Left and Right times

                        System.out.printf("Bar Z and TBar Calculation -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n", zBar, tBar, phiBar);
                        
                        lastComponent = component; // Avoid re-printing for the same component
                        leftHitFound = false;
                        rightHitFound = false;
                    }
                } else if (layer >= 10 && layer <= 19) { // Wedge hit
                    int wedgeIndex = component % 10; // Wedge index within the bar (0-9)
                    double zWedge = (wedgeIndex - (10 - 1) / 2.0) * WEDGE_SPACING; // Calculate ZWedge
                    double phiWedge = phiBar; // Same phi as the bar
                    double timeWedge = time;

                    System.out.printf("Wedge Hit -> Sector: %d, Layer: %d, Component: %d, ADC: %d, Time Wedge: %.2f ns, Pedestal: %d, Z Position: %.2f mm, Phi: %.2f rad\n",
                                      sector, layer, wedgeIndex, adc, time, pedestal, zWedge, phiWedge);

                    // Clustering
                    if (lastComponent == component) { // Check clustering conditions with the previous bar component
                        double deltaZ = Math.abs(zBar - zWedge);
                        double deltaPhi = Math.abs(phiBar - phiWedge);
                        double deltaTime = Math.abs(tBar - timeWedge);

                        if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                            clusterCount++;
                            int clusterSize = 3; // Minimum of 3 hits: 2 bar hits + 1 wedge hit
                            
                            // Updating the cluster size map
                            clusterSizeCount.put(clusterSize, clusterSizeCount.getOrDefault(clusterSize, 0) + 1);

                            System.out.println("Valid Cluster Formed:");
                            System.out.printf("  Bar Hits -> Sector: %d, Component: %d\n", sector, component);
                            System.out.printf("    Bar Hit -> ADC: %d, ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n",
                                              adc, zBar, tBar, phiBar);
                            System.out.printf("  Wedge Hit -> Sector: %d, Layer: %d, Component: %d, ADC: %d, Time Wedge: %.2f ns, Z Position: %.2f mm, Phi: %.2f rad\n",
                                              sector, layer, wedgeIndex, adc, timeWedge, zWedge, phiWedge);
                            System.out.printf("  Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", deltaZ, deltaPhi, deltaTime);
                        }
                    }
                }
            }
        }

        // Print cluster size distribution
        System.out.printf("\nTotal Valid Clusters Formed Across All Events: %d\n", clusterCount);
        for (Map.Entry<Integer, Integer> entry : clusterSizeCount.entrySet()) {
            System.out.printf("Clusters with %d hits: %d\n", entry.getKey(), entry.getValue());
        }
    }
}

*/


/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN{

    private static final double VEFF = 200.0; // mm/ns
    private static final double BAR_LENGTH = 280.0; // mm
    private static final double WEDGE_SPACING = 30.0; // mm
    private static final double Z_THRESHOLD = 150.0; // mm, threshold for clustering
    private static final double PHI_THRESHOLD = 0.1; // rad, threshold for clustering
    private static final double TIME_THRESHOLD = 1.0; // ns, threshold for clustering
    private static final int NUM_BARS = 60; // Number of bars
    private static final int NUM_WEDGES_PER_BAR = 10; // Wedges per bar

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            System.out.println("\nProcessing event with " + numHits + " hits...");

            boolean leftHitFound = false;
            boolean rightHitFound = false;
            double tBarLeft = 0.0;
            double tBarRight = 0.0;
            double zBar = 0.0;
            double tBar = 0.0;
            double phiBar = 0.0;
            int lastComponent = -1;

            for (int hitIndex = 0; hitIndex < numHits; hitIndex++) {
                int sector = atofAdcBank.getByte("sector", hitIndex);
                int layer = atofAdcBank.getByte("layer", hitIndex);
                int component = atofAdcBank.getShort("component", hitIndex);
                int order = atofAdcBank.getByte("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                float time = atofAdcBank.getFloat("time", hitIndex);
                int pedestal = atofAdcBank.getShort("ped", hitIndex);

                if (layer == 0) { // Bar hit
                    phiBar = 2 * Math.PI * component / NUM_BARS; // Phi calculation for the bar
                    if (order == 0) { // Left side
                        double zLeft = zBar - (BAR_LENGTH / 2);
                        tBarLeft = time - (zLeft / VEFF);
                        leftHitFound = true;
                        System.out.printf("Left Bar Hit -> Sector: %d, Layer: %d, Component: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, ZLeft: %.2f mm\n",
                                sector, layer, component, adc, time, pedestal, zLeft);
                    } else { // Right side
                        double zRight = zBar + (BAR_LENGTH / 2);
                        tBarRight = time - (zRight / VEFF);
                        rightHitFound = true;
                        System.out.printf("Right Bar Hit -> Sector: %d, Layer: %d, Component: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, ZRight: %.2f mm\n",
                                sector, layer, component, adc, time, pedestal, zRight);
                    }

                    if (leftHitFound && rightHitFound && lastComponent != component) { // Calculate ZBar and TBar after both hits
                        double deltaT = Math.abs(tBarLeft - tBarRight);
                        zBar = VEFF * deltaT / 2; // ZBar calculation
                        tBar = Math.min(tBarLeft, tBarRight); // Use the minimum of the times as TBar
                        System.out.printf("Bar Z and TBar Calculation -> ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n", zBar, tBar, phiBar);
                        
                        lastComponent = component; // Prevent re-processing the same component
                        leftHitFound = false;
                        rightHitFound = false;
                    }
                } else if (layer >= 10 && layer <= 19) { // Wedge hit
                    int wedgeIndex = component % NUM_WEDGES_PER_BAR; // Calculate wedge index (0-9)
                    double zWedge = (wedgeIndex - (NUM_WEDGES_PER_BAR - 1) / 2.0) * WEDGE_SPACING; // ZWedge calculation
                    double phiWedge = phiBar; // Same phi as the bar
                    double timeWedge = time;

                    System.out.printf("Wedge Hit -> Sector: %d, Layer: %d, Wedge Index: %d, ADC: %d, Time Wedge: %.2f ns, Pedestal: %d, Z Position: %.2f mm, Phi: %.2f rad\n",
                            sector, layer, wedgeIndex, adc, time, pedestal, zWedge, phiWedge);

                    // Clustering logic
                    if (lastComponent == component) {
                        double deltaZ = Math.abs(zBar - zWedge);
                        double deltaPhi = Math.abs(phiBar - phiWedge);
                        double deltaTime = Math.abs(tBar - timeWedge);

                        if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                            System.out.println("Valid Cluster Formed:");
                            System.out.printf("  Bar Hits -> Sector: %d, Component: %d\n", sector, component);
                            System.out.printf("    Bar Hit -> ADC (Left): %d, ADC (Right): %d, ZBar: %.2f mm, TBar: %.2f ns, Phi: %.2f rad\n",
                                    adc, adc, zBar, tBar, phiBar);
                            System.out.printf("  Wedge Hit -> Sector: %d, Layer: %d, Wedge Index: %d, ADC: %d, Time Wedge: %.2f ns, Z Position: %.2f mm, Phi: %.2f rad\n",
                                    sector, layer, wedgeIndex, adc, timeWedge, zWedge, phiWedge);
                            System.out.printf("  Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n", deltaZ, deltaPhi, deltaTime);
                        }
                    }
                }
            }
        }
    }
}

*/






/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

public class  ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final int NUM_WEDGES = 10;
    private static final int NUM_BARS = 60;
    private static final double WEDGE_SPACING = 30.0;
    private static final double VELOCITY_EFF = 200.0;

    private static final double Z_THRESHOLD = 5.0;
    private static final double TIME_THRESHOLD = 1.0;

    private static class EventData {
        Double zWedge = null, zBar = null, phiWedge = null, phiBar = null, timeWedge = null, timeBar = null;
        int sector, layer, component, order, adc;
        short pedestal;

        EventData(Double zWedge, Double zBar, Double phiWedge, Double phiBar,
                  Double timeWedge, Double timeBar, int sector, int layer,
                  int component, int order, int adc, short pedestal) {
            this.zWedge = zWedge;
            this.zBar = zBar;
            this.phiWedge = phiWedge;
            this.phiBar = phiBar;
            this.timeWedge = timeWedge;
            this.timeBar = timeBar;
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.pedestal = pedestal;
        }
    }

    private static class Cluster {
        List<EventData> events = new ArrayList<>();
        EventData wedgeHit = null;
        EventData barHitLeft = null;
        EventData barHitRight = null;

        public void addEvent(EventData event) {
            if (event.zBar != null && event.order == 0) {
                barHitLeft = event;
            } else if (event.zBar != null && event.order == 1) {
                barHitRight = event;
            } else if (event.zWedge != null) {
                wedgeHit = event;
            }
            events.add(event);
        }

        public boolean isValidCluster() {
            return (barHitLeft != null && barHitRight != null && wedgeHit != null);
        }

        public void printCluster(int clusterIndex) {
            System.out.printf("Cluster %d - Size: %d\n", clusterIndex, events.size());
            for (EventData event : events) {
                System.out.printf("  %s %s %s %s TW: %.2f, TB: %.2f, Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Pedestal: %d\n",
                        (event.zWedge != null ? String.format("ZW: %.2f", event.zWedge) : "ZW: N/A"),
                        (event.zBar != null ? String.format("ZB: %.2f", event.zBar) : "ZB: N/A"),
                        (event.phiWedge != null ? String.format("PhiW: %.2f", event.phiWedge) : "PhiW: N/A"),
                        (event.phiBar != null ? String.format("PhiB: %.2f", event.phiBar) : "PhiB: N/A"),
                        (event.timeWedge != null ? event.timeWedge : 0.0),
                        (event.timeBar != null ? event.timeBar : 0.0),
                        event.sector, event.layer, event.component, event.order, event.adc, event.pedestal);
            }

            // Calculation of delta values
            if (isValidCluster()) {
                double deltaZ = Math.abs(barHitLeft.zBar - wedgeHit.zWedge);
                double deltaPhi = Math.abs(barHitLeft.phiBar - wedgeHit.phiWedge);
                double deltaTime = Math.abs(barHitLeft.timeBar - wedgeHit.timeWedge);
                System.out.printf("Delta Z: %.2f, Delta Phi: %.2f, Delta Time: %.2f\n", deltaZ, deltaPhi, deltaTime);
            }
        }
    }

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();

        List<Integer> clusterSizes = new ArrayList<>();
        List<Integer> clusterIndices = new ArrayList<>();

        List<Double> deltaZList = new ArrayList<>();
        List<Double> deltaTimeList = new ArrayList<>();
        List<Double> deltaPhiList = new ArrayList<>();
        List<Integer> eventIndicesForDeltas = new ArrayList<>();

        int eventIndex = 0;
        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numRows = atofAdcBank.getRows();
            System.out.println("\nProcessing a new event with " + numRows + " hits...");

            List<EventData> eventsData = new ArrayList<>();
            List<Cluster> clusters = new ArrayList<>();

            Cluster currentCluster = new Cluster();

            for (int hitIndex = 0; hitIndex < numRows; hitIndex++) {
                int sector = atofAdcBank.getInt("sector", hitIndex);
                int layer = atofAdcBank.getInt("layer", hitIndex);
                int component = atofAdcBank.getInt("component", hitIndex);
                int order = atofAdcBank.getInt("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                short pedestal = atofAdcBank.getShort("ped", hitIndex);
                double time = atofAdcBank.getFloat("time", hitIndex);

                // Correct layer identification for wedge hits
                if (layer >= 16) {
                    layer = layer == 16 ? 2 : 3;
                }

                System.out.println("Hit -> Component: " + component + ", Layer: " + layer + ", Sector: " + sector + ", Order: " + order);

                if (component >= 1 && component <= NUM_BARS && (layer == 0 || layer == 1)) {
                    if (currentCluster.barHitLeft == null && order == 0) {
                        Double zBar = calculateZForBar(atofAdcBank, hitIndex);
                        double phiBar = calculatePhiForBar(component - 1);
                        currentCluster.addEvent(new EventData(null, zBar, null, phiBar, null, time, sector, layer, component, order, adc, pedestal));
                    } else if (currentCluster.barHitRight == null && order == 1) {
                        Double zBar = calculateZForBar(atofAdcBank, hitIndex);
                        double phiBar = calculatePhiForBar(component - 1);
                        currentCluster.addEvent(new EventData(null, zBar, null, phiBar, null, time, sector, layer, component, order, adc, pedestal));
                    }
                } else if (component >= 1 && component <= NUM_BARS && layer >= 2) {
                    int wedgeIndex = (component - 1) % NUM_WEDGES;
                    double wedgeZ = calculateZForWedge(wedgeIndex);
                    double wedgePhi = calculatePhiForBar(component - 1);
                    currentCluster.addEvent(new EventData(wedgeZ, null, wedgePhi, null, time, null, sector, layer, component, order, adc, pedestal));
                }
            }

            if (currentCluster.isValidCluster()) {
                clusters.add(currentCluster);
                System.out.println("Clusters formed in this event:");
                currentCluster.printCluster(clusters.size());

                deltaZList.add(currentCluster.barHitLeft.zBar - currentCluster.wedgeHit.zWedge);
                deltaPhiList.add(currentCluster.barHitLeft.phiBar - currentCluster.wedgeHit.phiWedge);
                deltaTimeList.add(currentCluster.barHitLeft.timeBar - currentCluster.wedgeHit.timeWedge);
                eventIndicesForDeltas.add(eventIndex);

                clusterSizes.add(currentCluster.events.size());
                clusterIndices.add(eventIndex);
            }

            eventIndex++;
        }

        plotClusterSizeVsEventIndex(clusterSizes, clusterIndices);
        plotDeltaVsEventIndex(deltaZList, eventIndicesForDeltas, "Delta Z vs Event Index", "Delta Z");
        plotDeltaVsEventIndex(deltaPhiList, eventIndicesForDeltas, "Delta Phi vs Event Index", "Delta Phi");
        plotDeltaVsEventIndex(deltaTimeList, eventIndicesForDeltas, "Delta Time vs Event Index", "Delta Time");
    }

    private static void plotClusterSizeVsEventIndex(List<Integer> clusterSizes, List<Integer> clusterIndices) {
        XYSeries series = new XYSeries("Cluster Size vs Event Index");
        for (int i = 0; i < clusterSizes.size(); i++) {
            series.add(clusterIndices.get(i), clusterSizes.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createScatterPlot(
                "Cluster Size vs Event Index",
                "Event Index", "Cluster Size",
                dataset, PlotOrientation.VERTICAL, true, true, false);

        displayChart(chart, "Cluster Size vs Event Index Plot");
    }

    private static void plotDeltaVsEventIndex(List<Double> deltaValues, List<Integer> eventIndices, String plotTitle, String yAxisLabel) {
        XYSeries series = new XYSeries(plotTitle);
        for (int i = 0; i < deltaValues.size(); i++) {
            series.add(eventIndices.get(i), deltaValues.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createScatterPlot(
                plotTitle,
                "Event Index", yAxisLabel,
                dataset, PlotOrientation.VERTICAL, true, true, false);

        displayChart(chart, plotTitle);
    }

    private static void displayChart(JFreeChart chart, String title) {
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(new ChartPanel(chart));
        frame.pack();
        frame.setVisible(true);
    }

    private static double calculateZForBar(Bank bank, int rowIndex) {
        double timeLeftPMT = bank.getFloat("time", rowIndex);
        double timeRightPMT = bank.getFloat("time", rowIndex + 1);

        double timeBar = Math.min(timeLeftPMT, timeRightPMT); // Use minimum time as timeBar
        double zBar = VELOCITY_EFF * (timeRightPMT - timeLeftPMT) / 2.0; // Calculate Z position based on time difference

        return timeBar; // Return timeBar as the minimum time
    }

    private static double calculateZForWedge(int wedgeIndex) {
        return (wedgeIndex - (NUM_WEDGES - 1) / 2.0) * WEDGE_SPACING;
    }

    private static double calculatePhiForBar(int barIndex) {
        double phiMin = -Math.PI;
        double phiMax = Math.PI;
        return phiMin + (phiMax - phiMin) * barIndex / NUM_BARS;
    }
}

*/




/*
package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final int NUM_WEDGES = 10;  // 10 wedges per bar
    private static final int NUM_BARS = 60;    // 60 bars total
    private static final double WEDGE_SPACING = 30.0;  // Wedge spacing in mm
    private static final double VELOCITY_EFF = 200.0;  // Effective velocity for time calculations

    private static final double Z_THRESHOLD = 5.0;
    private static final double TIME_THRESHOLD = 1.0;

    private static class EventData {
        Double zWedge = null, zBar = null, phiWedge = null, phiBar = null, timeWedge = null, timeBar = null;
        int sector, layer, component, order, adc;
        short pedestal;

        EventData(Double zWedge, Double zBar, Double phiWedge, Double phiBar,
                  Double timeWedge, Double timeBar, int sector, int layer,
                  int component, int order, int adc, short pedestal) {
            this.zWedge = zWedge;
            this.zBar = zBar;
            this.phiWedge = phiWedge;
            this.phiBar = phiBar;
            this.timeWedge = timeWedge;
            this.timeBar = timeBar;
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.pedestal = pedestal;
        }
    }

    private static class Cluster {
        List<EventData> events = new ArrayList<>();
        EventData wedgeHit = null;
        EventData barHitLeft = null;
        EventData barHitRight = null;

        public void addEvent(EventData event) {
            if (event.zBar != null && event.order == 0) {
                barHitLeft = event;
            } else if (event.zBar != null && event.order == 1) {
                barHitRight = event;
            } else if (event.zWedge != null) {
                wedgeHit = event;
            }
            events.add(event);
        }

        public boolean isValidCluster() {
            return (barHitLeft != null && barHitRight != null && wedgeHit != null);
        }

        public void printCluster(int clusterIndex) {
            System.out.printf("Cluster %d - Size: %d\n", clusterIndex, events.size());
            for (EventData event : events) {
                System.out.printf("  %s %s %s %s TW: %.2f, TB: %.2f, Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Pedestal: %d\n",
                        (event.zWedge != null ? String.format("ZW: %.2f", event.zWedge) : "ZW: N/A"),
                        (event.zBar != null ? String.format("ZB: %.2f", event.zBar) : "ZB: N/A"),
                        (event.phiWedge != null ? String.format("PhiW: %.2f", event.phiWedge) : "PhiW: N/A"),
                        (event.phiBar != null ? String.format("PhiB: %.2f", event.phiBar) : "PhiB: N/A"),
                        (event.timeWedge != null ? event.timeWedge : 0.0),
                        (event.timeBar != null ? event.timeBar : 0.0),
                        event.sector, event.layer, event.component, event.order, event.adc, event.pedestal);
            }

            // Calculate delta values
            if (isValidCluster()) {
                double deltaZ = Math.abs(barHitLeft.zBar - wedgeHit.zWedge);  // Z of the bar vs Z of the wedge
                double deltaPhi = Math.abs(barHitLeft.phiBar - wedgeHit.phiWedge);  // Phi of the bar vs Phi of the wedge
                double deltaTime = Math.abs(barHitLeft.timeBar - wedgeHit.timeWedge);  // Time of the bar vs Time of the wedge
                System.out.printf("Delta Z: %.2f, Delta Phi: %.2f, Delta Time: %.2f\n", deltaZ, deltaPhi, deltaTime);
            }
        }
    }

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        processEvents(reader);
        reader.close();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();

        List<Integer> clusterSizes = new ArrayList<>();
        List<Integer> clusterIndices = new ArrayList<>();

        List<Double> deltaZList = new ArrayList<>();
        List<Double> deltaTimeList = new ArrayList<>();
        List<Double> deltaPhiList = new ArrayList<>();
        List<Integer> eventIndicesForDeltas = new ArrayList<>();

        int eventIndex = 0;
        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numRows = atofAdcBank.getRows();
            System.out.println("\nProcessing a new event with " + numRows + " hits...");

            List<EventData> eventsData = new ArrayList<>();
            List<Cluster> clusters = new ArrayList<>();

            Cluster currentCluster = new Cluster();

            for (int hitIndex = 0; hitIndex < numRows; hitIndex++) {
                int sector = atofAdcBank.getInt("sector", hitIndex);
                int layer = atofAdcBank.getInt("layer", hitIndex);
                int component = atofAdcBank.getInt("component", hitIndex);
                int order = atofAdcBank.getInt("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                short pedestal = atofAdcBank.getShort("ped", hitIndex);
                double time = atofAdcBank.getFloat("time", hitIndex);

                // Correct layer identification for wedge hits (layers should be 2 or 3)
                if (layer >= 16) {
                    layer = layer == 16 ? 2 : 3;
                }

                System.out.println("Hit -> Component: " + component + ", Layer: " + layer + ", Sector: " + sector + ", Order: " + order);

                // Bar hit identification: Component range is 1 to 60, and layer 0 or 1 is for bars
                if (component >= 1 && component <= NUM_BARS && (layer == 0 || layer == 1)) {
                    if (currentCluster.barHitLeft == null && order == 0) {
                        Double zBar = calculateZForBar(atofAdcBank, hitIndex);
                        double phiBar = calculatePhiForBar(component - 1);
                        currentCluster.addEvent(new EventData(null, zBar, null, phiBar, null, time, sector, layer, component, order, adc, pedestal));
                    } else if (currentCluster.barHitRight == null && order == 1) {
                        Double zBar = calculateZForBar(atofAdcBank, hitIndex);
                        double phiBar = calculatePhiForBar(component - 1);
                        currentCluster.addEvent(new EventData(null, zBar, null, phiBar, null, time, sector, layer, component, order, adc, pedestal));
                    }
                } 
                // Wedge hit identification: Component range is 1 to 60 but layer 2 or 3 is for wedges
                else if (component >= 1 && component <= NUM_BARS && layer >= 2) {
                    int wedgeIndex = (component - 1) % NUM_WEDGES;
                    double wedgeZ = calculateZForWedge(wedgeIndex);
                    double wedgePhi = calculatePhiForBar(component - 1);
                    currentCluster.addEvent(new EventData(wedgeZ, null, wedgePhi, null, time, null, sector, layer, component, order, adc, pedestal));
                }
            }

            // Only form clusters if we have two bar hits and one wedge hit
            if (currentCluster.isValidCluster()) {
                clusters.add(currentCluster);
                System.out.println("Clusters formed in this event:");
                currentCluster.printCluster(clusters.size());

                // Collect delta values
                deltaZList.add(currentCluster.barHitLeft.zBar - currentCluster.wedgeHit.zWedge);
                deltaPhiList.add(currentCluster.barHitLeft.phiBar - currentCluster.wedgeHit.phiWedge);
                deltaTimeList.add(currentCluster.barHitLeft.timeBar - currentCluster.wedgeHit.timeWedge);
                eventIndicesForDeltas.add(eventIndex);

                // Add cluster size for plotting
                clusterSizes.add(currentCluster.events.size());
                clusterIndices.add(eventIndex);
            }

            eventIndex++;
        }

        // Plot cluster sizes
        plotClusterSizeVsEventIndex(clusterSizes, clusterIndices);

        // Plot delta Z, delta Phi, delta Time
        plotDeltaVsEventIndex(deltaZList, eventIndicesForDeltas, "Delta Z vs Event Index", "Delta Z");
        plotDeltaVsEventIndex(deltaPhiList, eventIndicesForDeltas, "Delta Phi vs Event Index", "Delta Phi");
        plotDeltaVsEventIndex(deltaTimeList, eventIndicesForDeltas, "Delta Time vs Event Index", "Delta Time");
    }

    private static void plotClusterSizeVsEventIndex(List<Integer> clusterSizes, List<Integer> clusterIndices) {
        XYSeries series = new XYSeries("Cluster Size vs Event Index");
        for (int i = 0; i < clusterSizes.size(); i++) {
            series.add(clusterIndices.get(i), clusterSizes.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createScatterPlot(
                "Cluster Size vs Event Index",
                "Event Index", "Cluster Size",
                dataset, PlotOrientation.VERTICAL, true, true, false);

        displayChart(chart, "Cluster Size vs Event Index Plot");
    }

    private static void plotDeltaVsEventIndex(List<Double> deltaValues, List<Integer> eventIndices, String plotTitle, String yAxisLabel) {
        XYSeries series = new XYSeries(plotTitle);
        for (int i = 0; i < deltaValues.size(); i++) {
            series.add(eventIndices.get(i), deltaValues.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createScatterPlot(
                plotTitle,
                "Event Index", yAxisLabel,
                dataset, PlotOrientation.VERTICAL, true, true, false);

        displayChart(chart, plotTitle);
    }

    private static void displayChart(JFreeChart chart, String title) {
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(new ChartPanel(chart));
        frame.pack();
        frame.setVisible(true);
    }

    private static double calculateZForBar(Bank bank, int rowIndex) {
        double timeLeftPMT = bank.getFloat("time", rowIndex);
        double timeRightPMT = bank.getFloat("time", rowIndex + 1);
        return VELOCITY_EFF * (timeRightPMT - timeLeftPMT) / 2.0;
    }

    private static double calculateZForWedge(int wedgeIndex) {
        return (wedgeIndex - (NUM_WEDGES - 1) / 2.0) * WEDGE_SPACING;
    }

    private static double calculatePhiForBar(int barIndex) {
        double phiMin = -Math.PI;
        double phiMax = Math.PI;
        return phiMin + (phiMax - phiMin) * barIndex / NUM_BARS;
    }
}
*/
