package edu.unc.mapseq.workflow.rnaseq;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.jgrapht.DirectedGraph;
import org.jgrapht.ext.VertexNameProvider;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.junit.Test;
import org.renci.jlrm.condor.ext.CondorDOTExporter;
import org.renci.jlrm.condor.CondorJob;
import org.renci.jlrm.condor.CondorJobEdge;

import edu.unc.mapseq.module.bedtools.CoverageBedCLI;
import edu.unc.mapseq.module.bowtie.BowtieBuildCLI;
import edu.unc.mapseq.module.core.CatCLI;
import edu.unc.mapseq.module.core.DetermineReadLengthCLI;
import edu.unc.mapseq.module.core.GUnZipCLI;
import edu.unc.mapseq.module.core.RegexCatCLI;
import edu.unc.mapseq.module.core.RemoveCLI;
import edu.unc.mapseq.module.core.SedCLI;
import edu.unc.mapseq.module.core.SortCLI;
import edu.unc.mapseq.module.filter.PruneISOFormsFromGeneQuantFileCLI;
import edu.unc.mapseq.module.filter.StripTrailingTabsCLI;
import edu.unc.mapseq.module.mapsplice.AlignmentHandlerMultiCLI;
import edu.unc.mapseq.module.mapsplice.ClusterCLI;
import edu.unc.mapseq.module.mapsplice.Filter1HitsCLI;
import edu.unc.mapseq.module.mapsplice.FilterJunctionByROCarguNonCanonicalCLI;
import edu.unc.mapseq.module.mapsplice.FilterOriginalFusionCLI;
import edu.unc.mapseq.module.mapsplice.FusionSAM2JunctionFilterAnchorNewFormatCLI;
import edu.unc.mapseq.module.mapsplice.JunctionDBFusionCLI;
import edu.unc.mapseq.module.mapsplice.JunctionSequenceConstructionCLI;
import edu.unc.mapseq.module.mapsplice.MapSpliceMultiThreadCLI;
import edu.unc.mapseq.module.mapsplice.NewSAM2JuntionCLI;
import edu.unc.mapseq.module.mapsplice.ParseClusterCLI;
import edu.unc.mapseq.module.mapsplice.RSEMCalculateExpressionCLI;
import edu.unc.mapseq.module.mapsplice.ReadChromoSizeCLI;
import edu.unc.mapseq.module.mapsplice.ReadsToUnmappedSAMCLI;
import edu.unc.mapseq.module.mapsplice.SetUnmappedBitFlagCLI;
import edu.unc.mapseq.module.picard.PicardAddOrReplaceReadGroupsCLI;
import edu.unc.mapseq.module.qc.NormBedExonQuantCLI;
import edu.unc.mapseq.module.qc.NormalizeQuartileCLI;
import edu.unc.mapseq.module.samtools.SAMToolsFlagstatCLI;
import edu.unc.mapseq.module.samtools.SAMToolsIndexCLI;
import edu.unc.mapseq.module.samtools.SAMToolsSortCLI;
import edu.unc.mapseq.module.samtools.SAMToolsViewCLI;
import edu.unc.mapseq.module.samtools.SortByReferenceAndNameCLI;
import edu.unc.mapseq.module.ubu.UBUFastqFormatterCLI;
import edu.unc.mapseq.module.ubu.UBUSamFilterCLI;
import edu.unc.mapseq.module.ubu.UBUSamJunctionCLI;
import edu.unc.mapseq.module.ubu.UBUSamTranslateCLI;

public class RNASeqWorkflowTest {

    @Test
    public void createDot() {

        DirectedGraph<CondorJob, CondorJobEdge> graph = new DefaultDirectedGraph<CondorJob, CondorJobEdge>(
                CondorJobEdge.class);

        int count = 0;

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
        CondorJob determineReadLengthJob = new CondorJob(String.format("%s_%d",
                DetermineReadLengthCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(determineReadLengthJob);
        graph.addEdge(fastqFormatterR1Job, determineReadLengthJob);
        graph.addEdge(fastqFormatterR2Job, determineReadLengthJob);

        // new job
        CondorJob readChromoSizesJob = new CondorJob(String.format("%s_%d", ReadChromoSizeCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(readChromoSizesJob);

        // new job
        CondorJob originalSAMMapspliceMultiThreadJob = new CondorJob(String.format("%s_%d",
                MapSpliceMultiThreadCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(originalSAMMapspliceMultiThreadJob);
        graph.addEdge(fastqFormatterR1Job, originalSAMMapspliceMultiThreadJob);
        graph.addEdge(fastqFormatterR2Job, originalSAMMapspliceMultiThreadJob);

        // new job
        CondorJob originalSAMRegexCatJob = new CondorJob(String.format("%s_%d", RegexCatCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(originalSAMRegexCatJob);
        graph.addEdge(originalSAMMapspliceMultiThreadJob, originalSAMRegexCatJob);

        // new job
        CondorJob sam2JunctionArrayJob = new CondorJob(String.format("%s_%d", NewSAM2JuntionCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(sam2JunctionArrayJob);
        graph.addEdge(determineReadLengthJob, sam2JunctionArrayJob);
        graph.addEdge(originalSAMRegexCatJob, sam2JunctionArrayJob);

        // new job
        CondorJob filter1HitsJob = new CondorJob(String.format("%s_%d", Filter1HitsCLI.class.getSimpleName(), ++count),
                null);
        graph.addVertex(filter1HitsJob);
        graph.addEdge(sam2JunctionArrayJob, filter1HitsJob);

        // new job
        CondorJob filterJunctionByROCarguNonCanonicalJob = new CondorJob(String.format("%s_%d",
                FilterJunctionByROCarguNonCanonicalCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(filterJunctionByROCarguNonCanonicalJob);
        graph.addEdge(filter1HitsJob, filterJunctionByROCarguNonCanonicalJob);

        // new job
        CondorJob junctionSequenceConstructionJob = new CondorJob(String.format("%s_%d",
                JunctionSequenceConstructionCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(junctionSequenceConstructionJob);
        graph.addEdge(filterJunctionByROCarguNonCanonicalJob, junctionSequenceConstructionJob);

        CondorJob bowtieBuildJob = new CondorJob(String.format("%s_%d", BowtieBuildCLI.class.getSimpleName(), ++count),
                null);
        graph.addVertex(bowtieBuildJob);
        graph.addEdge(junctionSequenceConstructionJob, bowtieBuildJob);

        // new job
        CondorJob remappedSAMMapspliceMultiThreadJob = new CondorJob(String.format("%s_%d",
                MapSpliceMultiThreadCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(remappedSAMMapspliceMultiThreadJob);
        graph.addEdge(bowtieBuildJob, remappedSAMMapspliceMultiThreadJob);

        // new job
        CondorJob remappedRegexCatJob = new CondorJob(
                String.format("%s_%d", RegexCatCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(remappedRegexCatJob);
        graph.addEdge(remappedSAMMapspliceMultiThreadJob, remappedRegexCatJob);

        // new job
        CondorJob remapUnmapped1RegexCatJob = new CondorJob(String.format("%s_%d", RegexCatCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(remapUnmapped1RegexCatJob);
        graph.addEdge(remappedSAMMapspliceMultiThreadJob, remapUnmapped1RegexCatJob);

        // new job
        CondorJob remapUnmapped2RegexCatJob = new CondorJob(String.format("%s_%d", RegexCatCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(remapUnmapped2RegexCatJob);
        graph.addEdge(remappedSAMMapspliceMultiThreadJob, remapUnmapped2RegexCatJob);

        // new job
        CondorJob alignmentHandlerMultiJob = new CondorJob(String.format("%s_%d",
                AlignmentHandlerMultiCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(alignmentHandlerMultiJob);
        graph.addEdge(determineReadLengthJob, alignmentHandlerMultiJob);
        graph.addEdge(remappedRegexCatJob, alignmentHandlerMultiJob);

        // new job
        CondorJob filteredNormalAlignmentsFusionPairedCatJob = new CondorJob(String.format("%s_%d",
                CatCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(filteredNormalAlignmentsFusionPairedCatJob);
        graph.addEdge(alignmentHandlerMultiJob, filteredNormalAlignmentsFusionPairedCatJob);

        // new job
        CondorJob filteredNormalAlignmentsSingleCatJob = new CondorJob(String.format("%s_%d",
                CatCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(filteredNormalAlignmentsSingleCatJob);
        graph.addEdge(alignmentHandlerMultiJob, filteredNormalAlignmentsSingleCatJob);

        // new job
        CondorJob filteredNormalAlignmentsPairedCatJob = new CondorJob(String.format("%s_%d",
                CatCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(filteredNormalAlignmentsPairedCatJob);
        graph.addEdge(alignmentHandlerMultiJob, filteredNormalAlignmentsPairedCatJob);

        // new job
        CondorJob sedJob = new CondorJob(String.format("%s_%d", SedCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(sedJob);
        graph.addEdge(filteredNormalAlignmentsFusionPairedCatJob, sedJob);

        // new job
        CondorJob parseClusterJob = new CondorJob(
                String.format("%s_%d", ParseClusterCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(parseClusterJob);
        graph.addEdge(sedJob, parseClusterJob);

        // new job
        CondorJob clusterJob = new CondorJob(String.format("%s_%d", ClusterCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(clusterJob);
        graph.addEdge(parseClusterJob, clusterJob);

        // new job
        CondorJob mapspliceMultiThreadFusionJob = new CondorJob(String.format("%s_%d",
                MapSpliceMultiThreadCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(mapspliceMultiThreadFusionJob);
        graph.addEdge(clusterJob, mapspliceMultiThreadFusionJob);
        graph.addEdge(remapUnmapped1RegexCatJob, mapspliceMultiThreadFusionJob);
        graph.addEdge(remapUnmapped2RegexCatJob, mapspliceMultiThreadFusionJob);

        // new job
        CondorJob fusionSAM2JunctionFilterAnchorNewFormatJob = new CondorJob(String.format("%s_%d",
                FusionSAM2JunctionFilterAnchorNewFormatCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(fusionSAM2JunctionFilterAnchorNewFormatJob);
        graph.addEdge(mapspliceMultiThreadFusionJob, fusionSAM2JunctionFilterAnchorNewFormatJob);

        // new job
        CondorJob filterOriginalFusionJob = new CondorJob(String.format("%s_%d",
                FilterOriginalFusionCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(filterOriginalFusionJob);
        graph.addEdge(fusionSAM2JunctionFilterAnchorNewFormatJob, filterOriginalFusionJob);

        // new job
        CondorJob junctionDBFusionJob = new CondorJob(String.format("%s_%d", JunctionDBFusionCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(junctionDBFusionJob);
        graph.addEdge(filterOriginalFusionJob, junctionDBFusionJob);

        // new job
        CondorJob synFusionIndexBowtieBuildJob = new CondorJob(String.format("%s_%d",
                BowtieBuildCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(synFusionIndexBowtieBuildJob);
        graph.addEdge(junctionDBFusionJob, synFusionIndexBowtieBuildJob);

        // new job
        mapspliceMultiThreadFusionJob = new CondorJob(String.format("%s_%d",
                MapSpliceMultiThreadCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(mapspliceMultiThreadFusionJob);
        graph.addEdge(synFusionIndexBowtieBuildJob, mapspliceMultiThreadFusionJob);

        // new job
        CondorJob sortJob = new CondorJob(String.format("%s_%d", SortCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(sortJob);
        graph.addEdge(mapspliceMultiThreadFusionJob, sortJob);

        // new job
        CondorJob fusionUnmapped1RegexCatJob = new CondorJob(String.format("%s_%d", RegexCatCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(fusionUnmapped1RegexCatJob);
        graph.addEdge(mapspliceMultiThreadFusionJob, fusionUnmapped1RegexCatJob);

        // new job
        CondorJob fusionUnmapped2RegexCatJob = new CondorJob(String.format("%s_%d", RegexCatCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(fusionUnmapped2RegexCatJob);
        graph.addEdge(mapspliceMultiThreadFusionJob, fusionUnmapped2RegexCatJob);

        // new job
        CondorJob fusionUnmapped1ReadsToUnmappedSAMJob = new CondorJob(String.format("%s_%d",
                ReadsToUnmappedSAMCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(fusionUnmapped1ReadsToUnmappedSAMJob);
        graph.addEdge(fusionUnmapped1RegexCatJob, fusionUnmapped1ReadsToUnmappedSAMJob);

        // new job
        CondorJob fusionUnmapped2ReadsToUnmappedSAMJob = new CondorJob(String.format("%s_%d",
                ReadsToUnmappedSAMCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(fusionUnmapped2ReadsToUnmappedSAMJob);
        graph.addEdge(fusionUnmapped2RegexCatJob, fusionUnmapped2ReadsToUnmappedSAMJob);

        // new job
        alignmentHandlerMultiJob = new CondorJob(String.format("%s_%d", AlignmentHandlerMultiCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(alignmentHandlerMultiJob);
        graph.addEdge(sortJob, alignmentHandlerMultiJob);
        graph.addEdge(determineReadLengthJob, alignmentHandlerMultiJob);

        // new job
        sortJob = new CondorJob(String.format("%s_%d", SortCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(sortJob);
        graph.addEdge(alignmentHandlerMultiJob, sortJob);
        graph.addEdge(fusionUnmapped1ReadsToUnmappedSAMJob, sortJob);
        graph.addEdge(fusionUnmapped2ReadsToUnmappedSAMJob, sortJob);

        // new job
        CondorJob setUnmappedBitFlagJob = new CondorJob(String.format("%s_%d",
                SetUnmappedBitFlagCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(setUnmappedBitFlagJob);
        graph.addEdge(sortJob, setUnmappedBitFlagJob);

        // new job
        CondorJob finalAlignmentsHeadedCatJob = new CondorJob(String.format("%s_%d", CatCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(finalAlignmentsHeadedCatJob);
        graph.addEdge(setUnmappedBitFlagJob, finalAlignmentsHeadedCatJob);
        graph.addEdge(readChromoSizesJob, finalAlignmentsHeadedCatJob);
        graph.addEdge(filteredNormalAlignmentsPairedCatJob, finalAlignmentsHeadedCatJob);
        graph.addEdge(filteredNormalAlignmentsSingleCatJob, finalAlignmentsHeadedCatJob);
        graph.addEdge(filteredNormalAlignmentsFusionPairedCatJob, finalAlignmentsHeadedCatJob);
        graph.addEdge(fusionUnmapped1ReadsToUnmappedSAMJob, finalAlignmentsHeadedCatJob);
        graph.addEdge(fusionUnmapped2ReadsToUnmappedSAMJob, finalAlignmentsHeadedCatJob);
        graph.addEdge(alignmentHandlerMultiJob, finalAlignmentsHeadedCatJob);

        // new job
        CondorJob samtoolsViewJob = new CondorJob(
                String.format("%s_%d", SAMToolsViewCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(samtoolsViewJob);
        graph.addEdge(finalAlignmentsHeadedCatJob, samtoolsViewJob);

        // new job
        CondorJob picardAddOrReplaceReadGroupsJob = new CondorJob(String.format("%s_%d",
                PicardAddOrReplaceReadGroupsCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(picardAddOrReplaceReadGroupsJob);
        graph.addEdge(samtoolsViewJob, picardAddOrReplaceReadGroupsJob);

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
        CondorJob ubuSamJunctionJob = new CondorJob(String.format("%s_%d", UBUSamJunctionCLI.class.getSimpleName(),
                ++count), null);

        graph.addVertex(ubuSamJunctionJob);
        graph.addEdge(samtoolsIndexJob, ubuSamJunctionJob);

        // new job
        CondorJob coverageBedJob = new CondorJob(String.format("%s_%d", CoverageBedCLI.class.getSimpleName(), ++count),
                null);
        graph.addVertex(coverageBedJob);
        graph.addEdge(samtoolsSortJob, coverageBedJob);

        // new job
        CondorJob normBedExonQuantJob = new CondorJob(String.format("%s_%d", NormBedExonQuantCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(normBedExonQuantJob);
        graph.addEdge(coverageBedJob, normBedExonQuantJob);

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

        // new job
        CondorJob removeJob = new CondorJob(String.format("%s_%d", RemoveCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(removeJob);
        graph.addEdge(normalizeGeneQuantJob, removeJob);
        graph.addEdge(normalizeISOFormQuantJob, removeJob);

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
        File dotFile = new File(srcSiteResourcesImagesDir, "workflow.dag.dot");
        try {
            FileWriter fw = new FileWriter(dotFile);
            dotExporter.export(fw, graph);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}
