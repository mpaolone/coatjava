package org.example;

public class Hit {
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

    public int getSector() {
        return sector;
    }

    public int getLayer() {
        return layer;
    }

    public int getComponent() {
        return component;
    }

    public int getOrder() {
        return order;
    }

    public int getAdc() {
        return adc;
    }

    public double getTime() {
        return time;
    }

    public double getPhi() {
        return phi;
    }

    public int getPedestal() {
        return pedestal;
    }
}


/*
package org.example;

public class Hit {
    private int sector;
    private int component;
    private double zPosition;
    private double time;
    private double phi;

    // Constructor
    public Hit(int sector, int component, double zPosition, double time, double phi) {
        this.sector = sector;
        this.component = component;
        this.zPosition = zPosition;
        this.time = time;
        this.phi = phi;
    }

    // Getter methods
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
*/
