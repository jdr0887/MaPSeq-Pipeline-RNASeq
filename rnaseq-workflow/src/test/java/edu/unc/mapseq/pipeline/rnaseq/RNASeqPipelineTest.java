package edu.unc.mapseq.pipeline.rnaseq;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.jgrapht.DirectedGraph;
import org.jgrapht.ext.VertexNameProvider;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.junit.Test;
import org.renci.jlrm.condor.CondorDOTExporter;
import org.renci.jlrm.condor.CondorJob;
import org.renci.jlrm.condor.CondorJobEdge;

import edu.unc.mapseq.module.core.GUnZipCLI;
import edu.unc.mapseq.module.core.MoveCLI;
import edu.unc.mapseq.module.filter.PruneISOFormsFromGeneQuantFileCLI;
import edu.unc.mapseq.module.filter.StripTrailingTabsCLI;
import edu.unc.mapseq.module.mapsplice.MapSpliceCLI;
import edu.unc.mapseq.module.mapsplice.RSEMCalculateExpressionCLI;
import edu.unc.mapseq.module.picard.PicardAddOrReplaceReadGroupsCLI;
import edu.unc.mapseq.module.qc.ExonQuantificationCLI;
import edu.unc.mapseq.module.qc.NormalizeQuartileCLI;
import edu.unc.mapseq.module.samtools.DetermineMedianCLI;
import edu.unc.mapseq.module.samtools.DetermineNumberOfReadsCLI;
import edu.unc.mapseq.module.samtools.SAMToolsFlagstatCLI;
import edu.unc.mapseq.module.samtools.SAMToolsIndexCLI;
import edu.unc.mapseq.module.samtools.SAMToolsPileupCLI;
import edu.unc.mapseq.module.samtools.SAMToolsSortCLI;
import edu.unc.mapseq.module.samtools.SortByReferenceAndNameCLI;
import edu.unc.mapseq.module.ubu.UBUFastqFormatterCLI;
import edu.unc.mapseq.module.ubu.UBUSamFilterCLI;
import edu.unc.mapseq.module.ubu.UBUSamJunctionCLI;
import edu.unc.mapseq.module.ubu.UBUSamTranslateCLI;

public class RNASeqPipelineTest {

    @Test
    public void createDot() {

        DirectedGraph<CondorJob, CondorJobEdge> graph = new DefaultDirectedGraph<CondorJob, CondorJobEdge>(
                CondorJobEdge.class);

        int count = 0;

        // new job
        CondorJob gunzipFastqR1Job = new CondorJob(String.format("%s_%d", GUnZipCLI.class.getSimpleName(), ++count),
                null);
        graph.addVertex(gunzipFastqR1Job);

        // new job
        CondorJob fastqFormatterR1Job = new CondorJob(String.format("%s_%d",
                UBUFastqFormatterCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(fastqFormatterR1Job);
        graph.addEdge(gunzipFastqR1Job, fastqFormatterR1Job);

        // new job
        CondorJob gunzipFastqR2Job = new CondorJob(String.format("%s_%d", GUnZipCLI.class.getSimpleName(), ++count),
                null);
        graph.addVertex(gunzipFastqR2Job);

        // new job
        CondorJob fastqFormatterR2Job = new CondorJob(String.format("%s_%d",
                UBUFastqFormatterCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(fastqFormatterR2Job);
        graph.addEdge(gunzipFastqR2Job, fastqFormatterR2Job);

        // new job
        CondorJob mapspliceJob = new CondorJob(String.format("%s_%d", MapSpliceCLI.class.getSimpleName(), ++count),
                null);
        graph.addVertex(mapspliceJob);
        graph.addEdge(fastqFormatterR1Job, mapspliceJob);
        graph.addEdge(fastqFormatterR2Job, mapspliceJob);

        // new job
        CondorJob moveJob = new CondorJob(String.format("%s_%d", MoveCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(moveJob);
        graph.addEdge(mapspliceJob, moveJob);

        // new job
        CondorJob picardAddOrReplaceReadGroupsJob = new CondorJob(String.format("%s_%d",
                PicardAddOrReplaceReadGroupsCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(picardAddOrReplaceReadGroupsJob);
        graph.addEdge(moveJob, picardAddOrReplaceReadGroupsJob);

        // new job
        CondorJob samtoolsSortJob = new CondorJob(
                String.format("%s_%d", SAMToolsSortCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(samtoolsSortJob);
        graph.addEdge(picardAddOrReplaceReadGroupsJob, samtoolsSortJob);

        // new job
        CondorJob samtoolsIndexJob = new CondorJob(String.format("%s_%d", SAMToolsIndexCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(samtoolsIndexJob);
        graph.addEdge(samtoolsSortJob, samtoolsIndexJob);

        // new job
        CondorJob samtoolsFlagstatJob = new CondorJob(String.format("%s_%d", SAMToolsFlagstatCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(samtoolsFlagstatJob);
        graph.addEdge(samtoolsIndexJob, samtoolsFlagstatJob);

        // new job
        CondorJob determineNumReadsJob = new CondorJob(String.format("%s_%d",
                DetermineNumberOfReadsCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(determineNumReadsJob);
        graph.addEdge(samtoolsFlagstatJob, determineNumReadsJob);

        // new job
        CondorJob ubuSamJunctionJob = new CondorJob(String.format("%s_%d", UBUSamJunctionCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(ubuSamJunctionJob);
        graph.addEdge(samtoolsIndexJob, ubuSamJunctionJob);

        // new job
        CondorJob determineMedianJob = new CondorJob(String.format("%s_%d", DetermineMedianCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(determineMedianJob);
        graph.addEdge(samtoolsSortJob, determineMedianJob);

        // new job
        CondorJob samtoolsPileupJob = new CondorJob(String.format("%s_%d", SAMToolsPileupCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(samtoolsPileupJob);
        graph.addEdge(samtoolsIndexJob, samtoolsPileupJob);

        // new job
        CondorJob exonQuantificationJob = new CondorJob(String.format("%s_%d",
                ExonQuantificationCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(exonQuantificationJob);
        graph.addEdge(samtoolsPileupJob, exonQuantificationJob);
        graph.addEdge(determineNumReadsJob, exonQuantificationJob);
        graph.addEdge(determineMedianJob, exonQuantificationJob);

        // new job
        CondorJob sortBAMByReferenceAndNameJob = new CondorJob(String.format("%s_%d",
                SortByReferenceAndNameCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(sortBAMByReferenceAndNameJob);
        graph.addEdge(samtoolsIndexJob, sortBAMByReferenceAndNameJob);

        // new job
        CondorJob ubuSamTranslateJob = new CondorJob(String.format("%s_%d", UBUSamTranslateCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(ubuSamTranslateJob);
        graph.addEdge(sortBAMByReferenceAndNameJob, ubuSamTranslateJob);

        // new job
        CondorJob ubuSamFilterJob = new CondorJob(
                String.format("%s_%d", UBUSamFilterCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(ubuSamFilterJob);
        graph.addEdge(ubuSamTranslateJob, ubuSamFilterJob);

        // new job
        CondorJob rsemJob = new CondorJob(String.format("%s_%d", RSEMCalculateExpressionCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(rsemJob);
        graph.addEdge(ubuSamFilterJob, rsemJob);

        // new job
        CondorJob rsemISOFormsResultStripTabJob = new CondorJob(String.format("%s_%d",
                StripTrailingTabsCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(rsemISOFormsResultStripTabJob);
        graph.addEdge(rsemJob, rsemISOFormsResultStripTabJob);

        // new job
        CondorJob pruneISOFormsFromGeneQuantFileJob = new CondorJob(String.format("%s_%d",
                PruneISOFormsFromGeneQuantFileCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(pruneISOFormsFromGeneQuantFileJob);
        graph.addEdge(rsemISOFormsResultStripTabJob, pruneISOFormsFromGeneQuantFileJob);

        // new job
        CondorJob normalizeGeneQuantJob = new CondorJob(String.format("%s_%d",
                NormalizeQuartileCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(normalizeGeneQuantJob);
        graph.addEdge(pruneISOFormsFromGeneQuantFileJob, normalizeGeneQuantJob);

        // new job
        CondorJob normalizeISOFormQuantJob = new CondorJob(String.format("%s_%d",
                NormalizeQuartileCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(normalizeISOFormQuantJob);
        graph.addEdge(pruneISOFormsFromGeneQuantFileJob, normalizeISOFormQuantJob);

        VertexNameProvider<CondorJob> vnpId = new VertexNameProvider<CondorJob>() {
            @Override
            public String getVertexName(CondorJob job) {
                return job.getName();
            }
        };

        VertexNameProvider<CondorJob> vnpLabel = new VertexNameProvider<CondorJob>() {
            @Override
            public String getVertexName(CondorJob job) {
                return job.getName();
            }
        };

        CondorDOTExporter<CondorJob, CondorJobEdge> dotExporter = new CondorDOTExporter<CondorJob, CondorJobEdge>(
                vnpId, vnpLabel, null, null, null, null);
        File srcSiteResourcesImagesDir = new File("src/site/resources/images");
        if (!srcSiteResourcesImagesDir.exists()) {
            srcSiteResourcesImagesDir.mkdirs();
        }
        File dotFile = new File(srcSiteResourcesImagesDir, "pipeline.dag.dot");
        try {
            FileWriter fw = new FileWriter(dotFile);
            dotExporter.export(fw, graph);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}
