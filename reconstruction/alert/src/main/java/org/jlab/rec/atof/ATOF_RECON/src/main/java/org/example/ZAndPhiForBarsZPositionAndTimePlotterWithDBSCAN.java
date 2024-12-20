
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

public class ZAndPhiForBarsZPositionAndTimePlotter{

    private static final int NUM_WEDGES = 10;  // Total number of wedges in each bar
    private static final int NUM_BARS = 60;    // Total number of bars
    private static final double WEDGE_SPACING = 30.0;  // Spacing between wedges in mm
    private static final double VELOCITY_EFF = 200.0;  // mm/ns, assumed effective velocity in the bar
    private static final double Z_THRESHOLD = 30.0;     // Maximum allowed difference in Z for clustering
    private static final double PHI_THRESHOLD = 0.01;  // Maximum allowed difference in Phi for clustering
    private static final double TIME_THRESHOLD = 1.0;  // Maximum allowed difference in Time for clustering

    // Class to represent event data for wedge and bar
    private static class EventData {
        double zWedge;
        double zBar;
        double phiWedge;
        double phiBar;
        double timeWedge;
        double timeBar;

        EventData(double zWedge, double zBar, double phiWedge, double phiBar, double timeWedge, double timeBar) {
            this.zWedge = zWedge;
            this.zBar = zBar;
            this.phiWedge = wrapPhi(phiWedge); // Ensure phi is within -π to π
            this.phiBar = wrapPhi(phiBar);     // Ensure phi is within -π to π
            this.timeWedge = timeWedge;
            this.timeBar = timeBar;
        }
    }

    // Class to represent clusters
    private static class Cluster {
        List<EventData> events = new ArrayList<>();

        public void addEvent(EventData event) {
            events.add(event);
        }

        public int getClusterSize() {
            return events.size();
        }

        public boolean isValidCluster() {
            return events.size() >= 2;  // Valid cluster must have at least 2 hits
        }

        public void printCluster(int clusterIndex) {
            System.out.printf("Cluster %d - Size: %d\n", clusterIndex, getClusterSize());
            for (EventData event : events) {
                System.out.printf("  ZW: %.2f, ZB: %.2f, PhiW: %.2f, PhiB: %.2f, TW: %.2f, TB: %.2f\n",
                        event.zWedge, event.zBar, event.phiWedge, event.phiBar, event.timeWedge, event.timeBar);
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

        List<EventData> eventsData = new ArrayList<>();
        List<Cluster> clusters = new ArrayList<>();

        // Extract and process wedge-bar association per event
        extractAndProcessEvents(reader, eventsData);

        // Form clusters based on proximity
        formClusters(eventsData, clusters);

        // Plot Delta Z, Delta Phi, Delta Time
        plotDeltas(eventsData);

        // Plot Cluster Size vs Event Pair Index or Cluster Index
        plotClusterSizeVsIndex(clusters);

        reader.close();
    }

    // Method to extract and process events per event and calculate deltas
    private static void extractAndProcessEvents(HipoReader reader, List<EventData> eventsData) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numRows = atofAdcBank.getRows();
            if (numRows < NUM_WEDGES) continue;  // Ensure there are enough rows for all wedges in the event

            System.out.println("\nProcessing a new event...");

            // Assume wedge and bar times and positions are retrieved or calculated as follows:
            for (int wedgeIndex = 0; wedgeIndex < NUM_WEDGES; wedgeIndex++) {
                // Calculate Z position for the wedge
                double wedgeZ = calculateZForWedge(wedgeIndex);

                // Get the Phi value based on the bar or sector (assuming equal spacing between bars)
                int barIndex = calculateBarIndex(atofAdcBank, wedgeIndex);
                double wedgePhi = calculatePhiForBar(barIndex);

                // Get wedge time from the ADC bank (assumed stored in the "time" column)
                double wedgeTime = atofAdcBank.getFloat("time", wedgeIndex);

                // Calculate Z, Phi, and time for the associated bar using the formula for Z
                double barZ = calculateZForBar(atofAdcBank, wedgeIndex);
                double barPhi = calculatePhiForBar(barIndex);
                double barTime = calculateBarTime(atofAdcBank, wedgeIndex); // Assuming geometric method for bar time

                // Store the wedge and bar data
                eventsData.add(new EventData(wedgeZ, barZ, wedgePhi, barPhi, wedgeTime, barTime));

                // Print the calculated Delta Z, Delta Phi, Delta Time
                double deltaZ = Math.abs(barZ - wedgeZ);
                double deltaPhi = Math.abs(barPhi - wedgePhi);
                double deltaTime = Math.abs(barTime - wedgeTime);
                System.out.printf("Event %d -> Delta Z: %.2f, Delta Phi: %.4f, Delta Time: %.4f\n", wedgeIndex, deltaZ, deltaPhi, deltaTime);
            }
        }
    }

    // **Updated** method to calculate Z for the bar based on time difference between left and right PMTs
    private static double calculateZForBar(Bank bank, int rowIndex) {
        double timeLeftPMT = bank.getFloat("time", rowIndex);     // Placeholder for left PMT time
        double timeRightPMT = bank.getFloat("time", rowIndex + 1); // Placeholder for right PMT time
        return VELOCITY_EFF * (timeRightPMT - timeLeftPMT) / 2.0;
    }

    // Method to calculate Z position for the wedge based on wedge index
    private static double calculateZForWedge(int wedgeIndex) {
        return (wedgeIndex - (NUM_WEDGES - 1) / 2.0) * WEDGE_SPACING;
    }

    // Calculate the Phi for a wedge (assuming even distribution of bars)
    private static double calculatePhiForBar(int barIndex) {
        double phi = -Math.PI + (2 * Math.PI) * barIndex / NUM_BARS;
        return wrapPhi(phi); // Ensure phi is within the range -π to π
    }

    // Ensure phi is between -π and π
    private static double wrapPhi(double phi) {
        while (phi > Math.PI) phi -= 2 * Math.PI;
        while (phi < -Math.PI) phi += 2 * Math.PI;
        return phi;
    }

    // Assuming a method to calculate bar time (e.g., using geometric method between two PMTs)
    private static double calculateBarTime(Bank bank, int rowIndex) {
        // Geometric method to get bar time using two PMT times (from ADC data)
        double timeLeftPMT = bank.getFloat("time", rowIndex);     // Placeholder for actual left PMT time
        double timeRightPMT = bank.getFloat("time", rowIndex + 1); // Placeholder for actual right PMT time
        return (timeLeftPMT + timeRightPMT) / 2;  // Midpoint of two PMTs
    }

    // Calculate bar index based on the sector, layer, and component (this logic may depend on your geometry)
    private static int calculateBarIndex(Bank bank, int rowIndex) {
        int sector = bank.getInt("sector", rowIndex);
        int layer = bank.getInt("layer", rowIndex);
        int component = bank.getInt("component", rowIndex);
        return sector * 4 * 60 + layer * 60 + component;  // Simple indexing logic
    }

    // Form clusters based on proximity in Z, Phi, and Time
    private static void formClusters(List<EventData> eventsData, List<Cluster> clusters) {
        for (EventData event : eventsData) {
            boolean addedToCluster = false;

            // Try to add the event to an existing cluster
            for (Cluster cluster : clusters) {
                if (isWithinProximity(cluster, event)) {
                    cluster.addEvent(event);
                    addedToCluster = true;
                    break;
                }
            }

            // If no suitable cluster was found, create a new one
            if (!addedToCluster) {
                Cluster newCluster = new Cluster();
                newCluster.addEvent(event);
                clusters.add(newCluster);
            }
        }

        // Print clusters
        int clusterIndex = 1;
        for (Cluster cluster : clusters) {
            if (cluster.isValidCluster()) {
                cluster.printCluster(clusterIndex++);
            }
        }
    }

    // Check if the event is within proximity of any event in the cluster
    private static boolean isWithinProximity(Cluster cluster, EventData event) {
        for (EventData clusterEvent : cluster.events) {
            if (Math.abs(clusterEvent.zBar - event.zBar) < Z_THRESHOLD &&
                    Math.abs(clusterEvent.phiBar - event.phiBar) < PHI_THRESHOLD &&
                    Math.abs(clusterEvent.timeBar - event.timeBar) < TIME_THRESHOLD) {
                return true;
            }
        }
        return false;
    }

    // Plot Delta Z, Delta Phi, and Delta Time
    private static void plotDeltas(List<EventData> eventsData) {
        XYSeries deltaZSeries = new XYSeries("Delta Z");
        XYSeries deltaPhiSeries = new XYSeries("Delta Phi");
        XYSeries deltaTimeSeries = new XYSeries("Delta Time");

        int eventIndex = 0; // Use event index for plotting

        for (EventData event : eventsData) {
            deltaZSeries.add(eventIndex, Math.abs(event.zBar - event.zWedge));
            deltaPhiSeries.add(eventIndex, Math.abs(event.phiBar - event.phiWedge));
            deltaTimeSeries.add(eventIndex, Math.abs(event.timeBar - event.timeWedge));
            eventIndex++;
        }

        plotSeries(deltaZSeries, "Delta Z", "Event Index", "Delta Z (mm)");
        plotSeries(deltaPhiSeries, "Delta Phi", "Event Index", "Delta Phi (radians)");
        plotSeries(deltaTimeSeries, "Delta Time", "Event Index", "Delta Time (ns)");
    }

    // Plot Cluster Size vs Event Pair Index or Cluster Index
    private static void plotClusterSizeVsIndex(List<Cluster> clusters) {
        XYSeries clusterSizeSeries = new XYSeries("Cluster Size");

        int clusterIndex = 1;
        for (Cluster cluster : clusters) {
            if (cluster.isValidCluster()) {
                System.out.println("Plotting cluster size: " + cluster.getClusterSize());
                clusterSizeSeries.add(clusterIndex++, cluster.getClusterSize());
            }
        }

        plotSeries(clusterSizeSeries, "Cluster Size", "Cluster Index", "Cluster Size (Number of Hits)");
    }

    // Utility method to display a plot for a given series
    private static void plotSeries(XYSeries series, String title, String xAxisLabel, String yAxisLabel) {
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset, PlotOrientation.VERTICAL, true, true, false);
        displayChart(chart, title);
    }

    // Utility method to display chart in a window
    private static void displayChart(JFreeChart chart, String title) {
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(new ChartPanel(chart));
        frame.pack();
        frame.setVisible(true);
    }
}
