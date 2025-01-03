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

public class ZAndPhiForBarsZPositionAndTimePlotter extends JFrame {

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

    public ZAndPhiForBarsZPositionAndTimePlotter(String title) {
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

      ZAndPhiForBarsZPositionAndTimePlotter demo = new ZAndPhiForBarsZPositionAndTimePlotter("Global Cluster Analysis Plots");
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



