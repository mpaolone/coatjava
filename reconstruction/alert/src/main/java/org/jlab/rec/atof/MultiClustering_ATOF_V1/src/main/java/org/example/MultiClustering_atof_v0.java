// This is the multi-clustering effort using ATOF::adc bank

package org.jlab.rec.atof.MultiClustering_ATOF_V0;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import java.util.ArrayList;
import java.util.List;

public class MultiClustering_atof_v0 {

    private static final double VEFF = 200.0;
    private static final double BAR_LENGTH = 280.0;
    private static final double WEDGE_SPACING = 30.0;
    private static final int N_WEDGE = 10;
    private static final double Z_THRESHOLD = 200.0;
    private static final double PHI_THRESHOLD = 0.2;
    private static final double TIME_THRESHOLD = 2.5;

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
        int eventId = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numHits = atofAdcBank.getRows();
            System.out.printf("\nProcessing Event ID: %d with %d hits...\n", eventId, numHits);

            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();
            extractHits(atofAdcBank, barHits, wedgeHits);

            List<Cluster> barClusters = formBarClusters(barHits);
            List<Cluster> barWedgeClusters = formBarWedgeClusters(barClusters, wedgeHits);

            int clusterId = 0;

            System.out.println("Bar-only Clusters:");
            for (Cluster cluster : barClusters) {
                System.out.printf("  Cluster ID: %d, Event ID: %d -> Z: %.2f mm, Phi: %.2f rad, Time: %.2f ns, Size: %d\n",
                        clusterId++, eventId, cluster.z, cluster.phi, cluster.time, cluster.hits.size());
                for (Hit hit : cluster.hits) {
                    System.out.printf("    Hit ID: %d -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Phi: %.2f rad\n",
                            hit.id, hit.sector, hit.layer, hit.component, hit.order, hit.adc, hit.time, hit.ped, hit.phi);
                }
            }

            System.out.println("Bar + Wedge Clusters:");
            for (Cluster cluster : barWedgeClusters) {
                System.out.printf("  Cluster ID: %d, Event ID: %d -> Z: %.2f mm, Phi: %.2f rad, Time: %.2f ns, Size: %d\n",
                        clusterId++, eventId, cluster.z, cluster.phi, cluster.time, cluster.hits.size());
                for (Hit hit : cluster.hits) {
                    String hitType = hit.layer == 0 ? "Bar" : "Wedge";
                    double zValue = hit.layer == 0 ? cluster.z : hit.z;
                    System.out.printf("    %s Hit ID: %d -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Time: %.2f ns, Pedestal: %d, Z: %.2f mm, Phi: %.2f rad\n",
                            hitType, hit.id, hit.sector, hit.layer, hit.component, hit.order, hit.adc, hit.time, hit.ped, zValue, hit.phi);
                }
                System.out.printf("    Thresholds -> Delta Z: %.2f mm, Delta Phi: %.2f rad, Delta Time: %.2f ns\n",
                        cluster.deltaZ, cluster.deltaPhi, cluster.deltaTime);
            }

            eventId++;
        }
    }

    private static void extractHits(Bank atofAdcBank, List<Hit> barHits, List<Hit> wedgeHits) {
        for (int i = 0; i < atofAdcBank.getRows(); i++) {
            int sector = atofAdcBank.getByte("sector", i);
            int layer = atofAdcBank.getByte("layer", i);
            int component = atofAdcBank.getShort("component", i);
            int order = atofAdcBank.getByte("order", i);
            int adc = atofAdcBank.getInt("ADC", i);
            float time = atofAdcBank.getFloat("time", i);
            int ped = atofAdcBank.getShort("ped", i);
            double phi = 2 * Math.PI * component / 60;
            double z = (layer == 0) ? 0.0 : (component % N_WEDGE - (N_WEDGE - 1) / 2.0) * WEDGE_SPACING;
            Hit hit = new Hit(sector, layer, component, order, adc, time, ped, z, phi);

            if (layer == 0) barHits.add(hit);
            else wedgeHits.add(hit);
        }
    }

    private static List<Cluster> formBarClusters(List<Hit> barHits) {
        List<Cluster> clusters = new ArrayList<>();
        for (int i = 0; i < barHits.size(); i++) {
            for (int j = i + 1; j < barHits.size(); j++) {
                Hit hit1 = barHits.get(i);
                Hit hit2 = barHits.get(j);
                if (hit1.component == hit2.component) {
                    double zBar = VEFF * (hit1.time - hit2.time) / 2;
                    double tMin = Math.min(hit1.time, hit2.time);
                    Cluster cluster = new Cluster(zBar, hit1.phi, tMin);
                    cluster.hits.add(hit1);
                    cluster.hits.add(hit2);
                    clusters.add(cluster);
                }
            }
        }
        return clusters;
    }

    private static List<Cluster> formBarWedgeClusters(List<Cluster> barClusters, List<Hit> wedgeHits) {
        List<Cluster> clusters = new ArrayList<>();
        for (Cluster barCluster : barClusters) {
            for (Hit wedgeHit : wedgeHits) {
                double deltaZ = Math.abs(barCluster.z - wedgeHit.z);
                double deltaPhi = Math.abs(barCluster.phi - wedgeHit.phi);
                double deltaTime = Math.abs(barCluster.time - wedgeHit.time);

                if (deltaZ < Z_THRESHOLD && deltaPhi < PHI_THRESHOLD && deltaTime < TIME_THRESHOLD) {
                    Cluster cluster = new Cluster(barCluster.z, barCluster.phi, Math.min(barCluster.time, wedgeHit.time));
                    cluster.hits.addAll(barCluster.hits);
                    cluster.hits.add(wedgeHit);
                    cluster.deltaZ = deltaZ;
                    cluster.deltaPhi = deltaPhi;
                    cluster.deltaTime = deltaTime;
                    clusters.add(cluster);
                }
            }
        }
        return clusters;
    }

    static class Hit {
        int sector, layer, component, order, adc, ped, id;
        double time, z, phi;

        Hit(int sector, int layer, int component, int order, int adc, double time, int ped, double z, double phi) {
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.time = time;
            this.ped = ped;
            this.z = z;
            this.phi = phi;
        }
    }

    static class Cluster {
        double z, phi, time, deltaZ, deltaPhi, deltaTime;
        List<Hit> hits = new ArrayList<>();

        Cluster(double z, double phi, double time) {
            this.z = z;
            this.phi = phi;
            this.time = time;
        }
    }
}

