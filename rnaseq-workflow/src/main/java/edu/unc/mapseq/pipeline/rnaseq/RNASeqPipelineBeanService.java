package edu.unc.mapseq.pipeline.rnaseq;

import edu.unc.mapseq.pipeline.AbstractPipelineBeanService;

public class RNASeqPipelineBeanService extends AbstractPipelineBeanService {

    private String bed;

    private String bowtieIndexDirectory;

    private String chromosomeDirectory;

    private String junctions;

    private String order;

    private String compositeExons;

    private String referenceSequencePrefix;

    private String siteName;

    public RNASeqPipelineBeanService() {
        super();
    }

    public String getSiteName() {
        return siteName;
    }

    public void setSiteName(String siteName) {
        this.siteName = siteName;
    }

    public String getBed() {
        return bed;
    }

    public void setBed(String bed) {
        this.bed = bed;
    }

    public String getBowtieIndexDirectory() {
        return bowtieIndexDirectory;
    }

    public void setBowtieIndexDirectory(String bowtieIndexDirectory) {
        this.bowtieIndexDirectory = bowtieIndexDirectory;
    }

    public String getChromosomeDirectory() {
        return chromosomeDirectory;
    }

    public void setChromosomeDirectory(String chromosomeDirectory) {
        this.chromosomeDirectory = chromosomeDirectory;
    }

    public String getJunctions() {
        return junctions;
    }

    public void setJunctions(String junctions) {
        this.junctions = junctions;
    }

    public String getOrder() {
        return order;
    }

    public void setOrder(String order) {
        this.order = order;
    }

    public String getCompositeExons() {
        return compositeExons;
    }

    public void setCompositeExons(String compositeExons) {
        this.compositeExons = compositeExons;
    }

    public String getReferenceSequencePrefix() {
        return referenceSequencePrefix;
    }

    public void setReferenceSequencePrefix(String referenceSequencePrefix) {
        this.referenceSequencePrefix = referenceSequencePrefix;
    }

}
