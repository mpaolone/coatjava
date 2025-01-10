package org.jlab.rec.atof.MultiClustering_ATOF_TDC_V2;

public class Hit {
    private int sector;
    private int component;
    private double zPosition;
    private double time;
    private double phi;

    
    public Hit(int sector, int component, double zPosition, double time, double phi) {
        this.sector = sector;
        this.component = component;
        this.zPosition = zPosition;
        this.time = time;
        this.phi = phi;
    }

    
    public int getSector() {
        return sector;
    }

    public int getComponent() {
        return component;
    }

    public double getzPosition() {
        return zPosition;
    }

    public double getTime() {
        return time;
    }

    public double getPhi() {
        return phi;
    }

    @Override
    public String toString() {
        return "Hit{Sector=" + sector + ", Component=" + component + ", Z=" + zPosition + 
               ", Time=" + time + ", Phi=" + phi + "}";
    }
}
