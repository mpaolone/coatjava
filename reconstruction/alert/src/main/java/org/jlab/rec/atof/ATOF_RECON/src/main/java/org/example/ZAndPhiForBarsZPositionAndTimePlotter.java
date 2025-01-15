package org.jlab.rec.atof.ATOF_RECON;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;

import java.util.ArrayList;
import java.util.List;

public class ZAndPhiForBarsZPositionAndTimePlotter {

    private static final int NUM_WEDGES = 10;  
    private static final int NUM_BARS = 60;    
    private static final double WEDGE_SPACING = 30.0; 
    private static final double VELOCITY_EFF = 200.0; 
    private static final double Z_THRESHOLD = 30.0;   
    private static final double PHI_THRESHOLD = 0.01; 
    private static final double TIME_THRESHOLD = 1.0; 

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
            this.phiWedge = wrapPhi(phiWedge);
            this.phiBar = wrapPhi(phiBar);
            this.timeWedge = timeWedge;
            this.timeBar = timeBar;
        }
    }

    private static class Cluster {
        List<EventData> events = new ArrayList<>();

        public void addEvent(EventData event) {
            events.add(event);
        }

        public int getClusterSize() {
            return events.size();
        }

        public boolean isValidCluster() {
            return events.size() >= 2;  
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

        extractAndProcessEvents(reader, eventsData);
        formClusters(eventsData, clusters);

        reader.close();
    }

    private static void extractAndProcessEvents(HipoReader reader, List<EventData> eventsData) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numRows = atofAdcBank.getRows();
            if (numRows < NUM_WEDGES) continue;

            System.out.println("\nProcessing a new event...");

            for (int wedgeIndex = 0; wedgeIndex < NUM_WEDGES; wedgeIndex++) {
                double wedgeZ = calculateZForWedge(wedgeIndex);
                int barIndex = calculateBarIndex(atofAdcBank, wedgeIndex);
                double wedgePhi = calculatePhiForBar(barIndex);
                double wedgeTime = atofAdcBank.getFloat("time", wedgeIndex);
                double barZ = calculateZForBar(atofAdcBank, wedgeIndex);
                double barPhi = calculatePhiForBar(barIndex);
                double barTime = calculateBarTime(atofAdcBank, wedgeIndex);

                eventsData.add(new EventData(wedgeZ, barZ, wedgePhi, barPhi, wedgeTime, barTime));

                double deltaZ = Math.abs(barZ - wedgeZ);
                double deltaPhi = Math.abs(barPhi - wedgePhi);
                double deltaTime = Math.abs(barTime - wedgeTime);
                System.out.printf("Event %d -> Delta Z: %.2f, Delta Phi: %.4f, Delta Time: %.4f\n", wedgeIndex, deltaZ, deltaPhi, deltaTime);
            }
        }
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
        double phi = -Math.PI + (2 * Math.PI) * barIndex / NUM_BARS;
        return wrapPhi(phi);
    }

    private static double wrapPhi(double phi) {
        while (phi > Math.PI) phi -= 2 * Math.PI;
        while (phi < -Math.PI) phi += 2 * Math.PI;
        return phi;
    }

    private static double calculateBarTime(Bank bank, int rowIndex) {
        double timeLeftPMT = bank.getFloat("time", rowIndex);
        double timeRightPMT = bank.getFloat("time", rowIndex + 1);
        return (timeLeftPMT + timeRightPMT) / 2;
    }

    private static int calculateBarIndex(Bank bank, int rowIndex) {
        int sector = bank.getInt("sector", rowIndex);
        int layer = bank.getInt("layer", rowIndex);
        int component = bank.getInt("component", rowIndex);
        return sector * 4 * 60 + layer * 60 + component;
    }

    private static void formClusters(List<EventData> eventsData, List<Cluster> clusters) {
        for (EventData event : eventsData) {
            boolean addedToCluster = false;

            for (Cluster cluster : clusters) {
                if (isWithinProximity(cluster, event)) {
                    cluster.addEvent(event);
                    addedToCluster = true;
                    break;
                }
            }

            if (!addedToCluster) {
                Cluster newCluster = new Cluster();
                newCluster.addEvent(event);
                clusters.add(newCluster);
            }
        }

        int clusterIndex = 1;
        for (Cluster cluster : clusters) {
            if (cluster.isValidCluster()) {
                cluster.printCluster(clusterIndex++);
            }
        }
    }

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
}

//TDC version of code

/*
package org.jlab.rec.atof.ATOF_RECON;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;

import java.util.ArrayList;
import java.util.List;

public class ZAndPhiForBarsZPositionAndTimePlotter_TDC {

    private static final int NUM_WEDGES = 10;
    private static final int NUM_BARS = 60;
    private static final double WEDGE_SPACING = 30.0;
    private static final double VELOCITY_EFF = 200.0;
    private static final double Z_THRESHOLD = 30.0;
    private static final double PHI_THRESHOLD = 0.01;
    private static final double TIME_THRESHOLD = 1.0;
    private static final double TDC_RESOLUTION = 0.015625; // TDC resolution in ns

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
            this.phiWedge = wrapPhi(phiWedge);
            this.phiBar = wrapPhi(phiBar);
            this.timeWedge = timeWedge;
            this.timeBar = timeBar;
        }
    }

    private static class Cluster {
        List<EventData> events = new ArrayList<>();

        public void addEvent(EventData event) {
            events.add(event);
        }

        public int getClusterSize() {
            return events.size();
        }

        public boolean isValidCluster() {
            return events.size() >= 2;
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

        if (!reader.getSchemaFactory().hasSchema("ATOF::tdc")) {
            System.err.println("Schema ATOF::tdc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        List<EventData> eventsData = new ArrayList<>();
        List<Cluster> clusters = new ArrayList<>();

        extractAndProcessEvents(reader, eventsData);
        formClusters(eventsData, clusters);

        reader.close();
    }

    private static void extractAndProcessEvents(HipoReader reader, List<EventData> eventsData) {
        Bank atofTdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::tdc"));
        Event event = new Event();

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofTdcBank);

            int numRows = atofTdcBank.getRows();
            if (numRows < NUM_WEDGES) continue;

            System.out.println("\nProcessing a new event...");

            for (int wedgeIndex = 0; wedgeIndex < NUM_WEDGES; wedgeIndex++) {
                double wedgeZ = calculateZForWedge(wedgeIndex);
                int barIndex = calculateBarIndex(atofTdcBank, wedgeIndex);
                double wedgePhi = calculatePhiForBar(barIndex);
                double wedgeTime = atofTdcBank.getInt("TDC", wedgeIndex) * TDC_RESOLUTION;
                double barZ = calculateZForBar(atofTdcBank, wedgeIndex);
                double barPhi = calculatePhiForBar(barIndex);
                double barTime = calculateBarTime(atofTdcBank, wedgeIndex);

                eventsData.add(new EventData(wedgeZ, barZ, wedgePhi, barPhi, wedgeTime, barTime));

                double deltaZ = Math.abs(barZ - wedgeZ);
                double deltaPhi = Math.abs(barPhi - wedgePhi);
                double deltaTime = Math.abs(barTime - wedgeTime);
                System.out.printf("Event %d -> Delta Z: %.2f, Delta Phi: %.4f, Delta Time: %.4f\n", wedgeIndex, deltaZ, deltaPhi, deltaTime);
            }
        }
    }

    private static double calculateZForBar(Bank bank, int rowIndex) {
        double timeLeftPMT = bank.getInt("TDC", rowIndex) * TDC_RESOLUTION;
        double timeRightPMT = bank.getInt("TDC", rowIndex + 1) * TDC_RESOLUTION;
        return VELOCITY_EFF * (timeRightPMT - timeLeftPMT) / 2.0;
    }

    private static double calculateZForWedge(int wedgeIndex) {
        return (wedgeIndex - (NUM_WEDGES - 1) / 2.0) * WEDGE_SPACING;
    }

    private static double calculatePhiForBar(int barIndex) {
        double phi = -Math.PI + (2 * Math.PI) * barIndex / NUM_BARS;
        return wrapPhi(phi);
    }

    private static double wrapPhi(double phi) {
        while (phi > Math.PI) phi -= 2 * Math.PI;
        while (phi < -Math.PI) phi += 2 * Math.PI;
        return phi;
    }

    private static double calculateBarTime(Bank bank, int rowIndex) {
        double timeLeftPMT = bank.getInt("TDC", rowIndex) * TDC_RESOLUTION;
        double timeRightPMT = bank.getInt("TDC", rowIndex + 1) * TDC_RESOLUTION;
        return (timeLeftPMT + timeRightPMT) / 2;
    }

    private static int calculateBarIndex(Bank bank, int rowIndex) {
        int sector = bank.getInt("sector", rowIndex);
        int layer = bank.getInt("layer", rowIndex);
        int component = bank.getInt("component", rowIndex);
        return sector * 4 * 60 + layer * 60 + component;
    }

    private static void formClusters(List<EventData> eventsData, List<Cluster> clusters) {
        for (EventData event : eventsData) {
            boolean addedToCluster = false;

            for (Cluster cluster : clusters) {
                if (isWithinProximity(cluster, event)) {
                    cluster.addEvent(event);
                    addedToCluster = true;
                    break;
                }
            }

            if (!addedToCluster) {
                Cluster newCluster = new Cluster();
                newCluster.addEvent(event);
                clusters.add(newCluster);
            }
        }

        int clusterIndex = 1;
        for (Cluster cluster : clusters) {
            if (cluster.isValidCluster()) {
                cluster.printCluster(clusterIndex++);
            }
        }
    }

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
}
*/






