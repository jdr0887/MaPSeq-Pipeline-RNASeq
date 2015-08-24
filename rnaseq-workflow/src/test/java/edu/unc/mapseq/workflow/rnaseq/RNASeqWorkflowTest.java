package edu.unc.mapseq.workflow.rnaseq;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.jgrapht.DirectedGraph;
import org.jgrapht.ext.VertexNameProvider;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.junit.Test;
import org.renci.jlrm.condor.CondorJob;
import org.renci.jlrm.condor.CondorJobEdge;
import org.renci.jlrm.condor.ext.CondorDOTExporter;

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
import edu.unc.mapseq.workflow.WorkflowException;
import edu.unc.mapseq.workflow.impl.WorkflowJobFactory;

public class RNASeqWorkflowTest {

    @Test
    public void createDot() {

        DirectedGraph<CondorJob, CondorJobEdge> graph = new DefaultDirectedGraph<CondorJob, CondorJobEdge>(
                CondorJobEdge.class);

        int count = 0;

        try {
            CondorJob gunzipFastqR1Job = WorkflowJobFactory.createJob(++count, GUnZipCLI.class, null).build();
            graph.addVertex(gunzipFastqR1Job);

            // new job
            CondorJob fastqFormatterR1Job = WorkflowJobFactory.createJob(++count, UBUFastqFormatterCLI.class, null)
                    .build();
            graph.addVertex(fastqFormatterR1Job);
            graph.addEdge(gunzipFastqR1Job, fastqFormatterR1Job);

            // new job
            CondorJob gunzipFastqR2Job = WorkflowJobFactory.createJob(++count, GUnZipCLI.class, null).build();
            graph.addVertex(gunzipFastqR2Job);

            // new job
            CondorJob fastqFormatterR2Job = WorkflowJobFactory.createJob(++count, UBUFastqFormatterCLI.class, null)
                    .build();
            graph.addVertex(fastqFormatterR2Job);
            graph.addEdge(gunzipFastqR2Job, fastqFormatterR2Job);

            // new job
            CondorJob determineReadLengthJob = WorkflowJobFactory
                    .createJob(++count, DetermineReadLengthCLI.class, null).build();
            graph.addVertex(determineReadLengthJob);
            graph.addEdge(fastqFormatterR1Job, determineReadLengthJob);
            graph.addEdge(fastqFormatterR2Job, determineReadLengthJob);

            // new job
            CondorJob readChromoSizesJob = WorkflowJobFactory.createJob(++count, ReadChromoSizeCLI.class, null).build();
            graph.addVertex(readChromoSizesJob);

            // new job
            CondorJob originalSAMMapspliceMultiThreadJob = WorkflowJobFactory.createJob(++count,
                    MapSpliceMultiThreadCLI.class, null).build();
            graph.addVertex(originalSAMMapspliceMultiThreadJob);
            graph.addEdge(fastqFormatterR1Job, originalSAMMapspliceMultiThreadJob);
            graph.addEdge(fastqFormatterR2Job, originalSAMMapspliceMultiThreadJob);

            // new job
            CondorJob originalSAMRegexCatJob = WorkflowJobFactory.createJob(++count, RegexCatCLI.class, null).build();
            graph.addVertex(originalSAMRegexCatJob);
            graph.addEdge(originalSAMMapspliceMultiThreadJob, originalSAMRegexCatJob);

            // new job
            CondorJob sam2JunctionArrayJob = WorkflowJobFactory.createJob(++count, NewSAM2JuntionCLI.class, null)
                    .build();
            graph.addVertex(sam2JunctionArrayJob);
            graph.addEdge(determineReadLengthJob, sam2JunctionArrayJob);
            graph.addEdge(originalSAMRegexCatJob, sam2JunctionArrayJob);

            // new job
            CondorJob filter1HitsJob = WorkflowJobFactory.createJob(++count, Filter1HitsCLI.class, null).build();
            graph.addVertex(filter1HitsJob);
            graph.addEdge(sam2JunctionArrayJob, filter1HitsJob);

            // new job
            CondorJob filterJunctionByROCarguNonCanonicalJob = WorkflowJobFactory.createJob(++count,
                    FilterJunctionByROCarguNonCanonicalCLI.class, null).build();
            graph.addVertex(filterJunctionByROCarguNonCanonicalJob);
            graph.addEdge(filter1HitsJob, filterJunctionByROCarguNonCanonicalJob);

            // new job
            CondorJob junctionSequenceConstructionJob = WorkflowJobFactory.createJob(++count,
                    JunctionSequenceConstructionCLI.class, null).build();
            graph.addVertex(junctionSequenceConstructionJob);
            graph.addEdge(filterJunctionByROCarguNonCanonicalJob, junctionSequenceConstructionJob);

            CondorJob bowtieBuildJob = WorkflowJobFactory.createJob(++count, BowtieBuildCLI.class, null).build();
            graph.addVertex(bowtieBuildJob);
            graph.addEdge(junctionSequenceConstructionJob, bowtieBuildJob);

            // new job
            CondorJob remappedSAMMapspliceMultiThreadJob = WorkflowJobFactory.createJob(++count,
                    MapSpliceMultiThreadCLI.class, null).build();
            graph.addVertex(remappedSAMMapspliceMultiThreadJob);
            graph.addEdge(bowtieBuildJob, remappedSAMMapspliceMultiThreadJob);

            // new job
            CondorJob remappedRegexCatJob = WorkflowJobFactory.createJob(++count, RegexCatCLI.class, null).build();
            graph.addVertex(remappedRegexCatJob);
            graph.addEdge(remappedSAMMapspliceMultiThreadJob, remappedRegexCatJob);

            // new job
            CondorJob remapUnmapped1RegexCatJob = WorkflowJobFactory.createJob(++count, RegexCatCLI.class, null)
                    .build();
            graph.addVertex(remapUnmapped1RegexCatJob);
            graph.addEdge(remappedSAMMapspliceMultiThreadJob, remapUnmapped1RegexCatJob);

            // new job
            CondorJob remapUnmapped2RegexCatJob = WorkflowJobFactory.createJob(++count, RegexCatCLI.class, null)
                    .build();
            graph.addVertex(remapUnmapped2RegexCatJob);
            graph.addEdge(remappedSAMMapspliceMultiThreadJob, remapUnmapped2RegexCatJob);

            // new job
            CondorJob alignmentHandlerMultiJob = WorkflowJobFactory.createJob(++count, AlignmentHandlerMultiCLI.class,
                    null).build();
            graph.addVertex(alignmentHandlerMultiJob);
            graph.addEdge(determineReadLengthJob, alignmentHandlerMultiJob);
            graph.addEdge(remappedRegexCatJob, alignmentHandlerMultiJob);

            // new job
            CondorJob filteredNormalAlignmentsFusionPairedCatJob = WorkflowJobFactory.createJob(++count, CatCLI.class,
                    null).build();
            graph.addVertex(filteredNormalAlignmentsFusionPairedCatJob);
            graph.addEdge(alignmentHandlerMultiJob, filteredNormalAlignmentsFusionPairedCatJob);

            // new job
            CondorJob filteredNormalAlignmentsSingleCatJob = WorkflowJobFactory.createJob(++count, CatCLI.class, null)
                    .build();
            graph.addVertex(filteredNormalAlignmentsSingleCatJob);
            graph.addEdge(alignmentHandlerMultiJob, filteredNormalAlignmentsSingleCatJob);

            // new job
            CondorJob filteredNormalAlignmentsPairedCatJob = WorkflowJobFactory.createJob(++count, CatCLI.class, null)
                    .build();
            graph.addVertex(filteredNormalAlignmentsPairedCatJob);
            graph.addEdge(alignmentHandlerMultiJob, filteredNormalAlignmentsPairedCatJob);

            // new job
            CondorJob sedJob = WorkflowJobFactory.createJob(++count, SedCLI.class, null).build();
            graph.addVertex(sedJob);
            graph.addEdge(filteredNormalAlignmentsFusionPairedCatJob, sedJob);

            // new job
            CondorJob parseClusterJob = WorkflowJobFactory.createJob(++count, ParseClusterCLI.class, null).build();
            graph.addVertex(parseClusterJob);
            graph.addEdge(sedJob, parseClusterJob);

            // new job
            CondorJob clusterJob = WorkflowJobFactory.createJob(++count, ClusterCLI.class, null).build();
            graph.addVertex(clusterJob);
            graph.addEdge(parseClusterJob, clusterJob);

            // new job
            CondorJob mapspliceMultiThreadFusionJob = WorkflowJobFactory.createJob(++count,
                    MapSpliceMultiThreadCLI.class, null).build();
            graph.addVertex(mapspliceMultiThreadFusionJob);
            graph.addEdge(clusterJob, mapspliceMultiThreadFusionJob);
            graph.addEdge(remapUnmapped1RegexCatJob, mapspliceMultiThreadFusionJob);
            graph.addEdge(remapUnmapped2RegexCatJob, mapspliceMultiThreadFusionJob);

            // new job
            CondorJob fusionSAM2JunctionFilterAnchorNewFormatJob = WorkflowJobFactory.createJob(++count,
                    FusionSAM2JunctionFilterAnchorNewFormatCLI.class, null).build();
            graph.addVertex(fusionSAM2JunctionFilterAnchorNewFormatJob);
            graph.addEdge(mapspliceMultiThreadFusionJob, fusionSAM2JunctionFilterAnchorNewFormatJob);

            // new job
            CondorJob filterOriginalFusionJob = WorkflowJobFactory.createJob(++count, FilterOriginalFusionCLI.class,
                    null).build();
            graph.addVertex(filterOriginalFusionJob);
            graph.addEdge(fusionSAM2JunctionFilterAnchorNewFormatJob, filterOriginalFusionJob);

            // new job
            CondorJob junctionDBFusionJob = WorkflowJobFactory.createJob(++count, JunctionDBFusionCLI.class, null)
                    .build();
            graph.addVertex(junctionDBFusionJob);
            graph.addEdge(filterOriginalFusionJob, junctionDBFusionJob);

            // new job
            CondorJob synFusionIndexBowtieBuildJob = WorkflowJobFactory.createJob(++count, BowtieBuildCLI.class, null)
                    .build();
            graph.addVertex(synFusionIndexBowtieBuildJob);
            graph.addEdge(junctionDBFusionJob, synFusionIndexBowtieBuildJob);

            // new job
            mapspliceMultiThreadFusionJob = WorkflowJobFactory.createJob(++count, MapSpliceMultiThreadCLI.class, null)
                    .build();
            graph.addVertex(mapspliceMultiThreadFusionJob);
            graph.addEdge(synFusionIndexBowtieBuildJob, mapspliceMultiThreadFusionJob);

            // new job
            CondorJob sortJob = WorkflowJobFactory.createJob(++count, SortCLI.class, null).build();
            graph.addVertex(sortJob);
            graph.addEdge(mapspliceMultiThreadFusionJob, sortJob);

            // new job
            CondorJob fusionUnmapped1RegexCatJob = WorkflowJobFactory.createJob(++count, RegexCatCLI.class, null)
                    .build();
            graph.addVertex(fusionUnmapped1RegexCatJob);
            graph.addEdge(mapspliceMultiThreadFusionJob, fusionUnmapped1RegexCatJob);

            // new job
            CondorJob fusionUnmapped2RegexCatJob = WorkflowJobFactory.createJob(++count, RegexCatCLI.class, null)
                    .build();
            graph.addVertex(fusionUnmapped2RegexCatJob);
            graph.addEdge(mapspliceMultiThreadFusionJob, fusionUnmapped2RegexCatJob);

            // new job
            CondorJob fusionUnmapped1ReadsToUnmappedSAMJob = WorkflowJobFactory.createJob(++count,
                    ReadsToUnmappedSAMCLI.class, null).build();
            graph.addVertex(fusionUnmapped1ReadsToUnmappedSAMJob);
            graph.addEdge(fusionUnmapped1RegexCatJob, fusionUnmapped1ReadsToUnmappedSAMJob);

            // new job
            CondorJob fusionUnmapped2ReadsToUnmappedSAMJob = WorkflowJobFactory.createJob(++count,
                    ReadsToUnmappedSAMCLI.class, null).build();
            graph.addVertex(fusionUnmapped2ReadsToUnmappedSAMJob);
            graph.addEdge(fusionUnmapped2RegexCatJob, fusionUnmapped2ReadsToUnmappedSAMJob);

            // new job
            alignmentHandlerMultiJob = WorkflowJobFactory.createJob(++count, AlignmentHandlerMultiCLI.class, null)
                    .build();
            graph.addVertex(alignmentHandlerMultiJob);
            graph.addEdge(sortJob, alignmentHandlerMultiJob);
            graph.addEdge(determineReadLengthJob, alignmentHandlerMultiJob);

            // new job
            sortJob = WorkflowJobFactory.createJob(++count, SortCLI.class, null).build();
            graph.addVertex(sortJob);
            graph.addEdge(alignmentHandlerMultiJob, sortJob);
            graph.addEdge(fusionUnmapped1ReadsToUnmappedSAMJob, sortJob);
            graph.addEdge(fusionUnmapped2ReadsToUnmappedSAMJob, sortJob);

            // new job
            CondorJob setUnmappedBitFlagJob = WorkflowJobFactory.createJob(++count, SetUnmappedBitFlagCLI.class, null)
                    .build();
            graph.addVertex(setUnmappedBitFlagJob);
            graph.addEdge(sortJob, setUnmappedBitFlagJob);

            // new job
            CondorJob finalAlignmentsHeadedCatJob = WorkflowJobFactory.createJob(++count, CatCLI.class, null).build();
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
            CondorJob samtoolsViewJob = WorkflowJobFactory.createJob(++count, SAMToolsViewCLI.class, null).build();
            graph.addVertex(samtoolsViewJob);
            graph.addEdge(finalAlignmentsHeadedCatJob, samtoolsViewJob);

            // new job
            CondorJob picardAddOrReplaceReadGroupsJob = WorkflowJobFactory.createJob(++count,
                    PicardAddOrReplaceReadGroupsCLI.class, null).build();
            graph.addVertex(picardAddOrReplaceReadGroupsJob);
            graph.addEdge(samtoolsViewJob, picardAddOrReplaceReadGroupsJob);

            // new job
            CondorJob samtoolsSortJob = WorkflowJobFactory.createJob(++count, SAMToolsSortCLI.class, null).build();
            graph.addVertex(samtoolsSortJob);
            graph.addEdge(picardAddOrReplaceReadGroupsJob, samtoolsSortJob);

            // new job
            CondorJob samtoolsIndexJob = WorkflowJobFactory.createJob(++count, SAMToolsIndexCLI.class, null).build();
            graph.addVertex(samtoolsIndexJob);
            graph.addEdge(samtoolsSortJob, samtoolsIndexJob);

            // new job
            CondorJob samtoolsFlagstatJob = WorkflowJobFactory.createJob(++count, SAMToolsFlagstatCLI.class, null)
                    .build();
            graph.addVertex(samtoolsFlagstatJob);
            graph.addEdge(samtoolsIndexJob, samtoolsFlagstatJob);

            // new job
            CondorJob ubuSamJunctionJob = WorkflowJobFactory.createJob(++count, UBUSamJunctionCLI.class, null).build();

            graph.addVertex(ubuSamJunctionJob);
            graph.addEdge(samtoolsIndexJob, ubuSamJunctionJob);

            // new job
            CondorJob coverageBedJob = WorkflowJobFactory.createJob(++count, CoverageBedCLI.class, null).build();
            graph.addVertex(coverageBedJob);
            graph.addEdge(samtoolsSortJob, coverageBedJob);

            // new job
            CondorJob normBedExonQuantJob = WorkflowJobFactory.createJob(++count, NormBedExonQuantCLI.class, null)
                    .build();
            graph.addVertex(normBedExonQuantJob);
            graph.addEdge(coverageBedJob, normBedExonQuantJob);

            // new job
            CondorJob sortBAMByReferenceAndNameJob = WorkflowJobFactory.createJob(++count,
                    SortByReferenceAndNameCLI.class, null).build();
            graph.addVertex(sortBAMByReferenceAndNameJob);
            graph.addEdge(samtoolsIndexJob, sortBAMByReferenceAndNameJob);

            // new job
            CondorJob ubuSamTranslateJob = WorkflowJobFactory.createJob(++count, UBUSamTranslateCLI.class, null)
                    .build();
            graph.addVertex(ubuSamTranslateJob);
            graph.addEdge(sortBAMByReferenceAndNameJob, ubuSamTranslateJob);

            // new job
            CondorJob ubuSamFilterJob = WorkflowJobFactory.createJob(++count, UBUSamFilterCLI.class, null).build();
            graph.addVertex(ubuSamFilterJob);
            graph.addEdge(ubuSamTranslateJob, ubuSamFilterJob);

            // new job
            CondorJob rsemJob = WorkflowJobFactory.createJob(++count, RSEMCalculateExpressionCLI.class, null).build();
            graph.addVertex(rsemJob);
            graph.addEdge(ubuSamFilterJob, rsemJob);

            // new job
            CondorJob rsemISOFormsResultStripTabJob = WorkflowJobFactory.createJob(++count, StripTrailingTabsCLI.class,
                    null).build();
            graph.addVertex(rsemISOFormsResultStripTabJob);
            graph.addEdge(rsemJob, rsemISOFormsResultStripTabJob);

            // new job
            CondorJob pruneISOFormsFromGeneQuantFileJob = WorkflowJobFactory.createJob(++count,
                    PruneISOFormsFromGeneQuantFileCLI.class, null).build();
            graph.addVertex(pruneISOFormsFromGeneQuantFileJob);
            graph.addEdge(rsemISOFormsResultStripTabJob, pruneISOFormsFromGeneQuantFileJob);

            // new job
            CondorJob normalizeGeneQuantJob = WorkflowJobFactory.createJob(++count, NormalizeQuartileCLI.class, null)
                    .build();
            graph.addVertex(normalizeGeneQuantJob);
            graph.addEdge(pruneISOFormsFromGeneQuantFileJob, normalizeGeneQuantJob);

            // new job
            CondorJob normalizeISOFormQuantJob = WorkflowJobFactory
                    .createJob(++count, NormalizeQuartileCLI.class, null).build();
            graph.addVertex(normalizeISOFormQuantJob);
            graph.addEdge(pruneISOFormsFromGeneQuantFileJob, normalizeISOFormQuantJob);

            // new job
            CondorJob removeJob = WorkflowJobFactory.createJob(++count, RemoveCLI.class, null).build();
            graph.addVertex(removeJob);
            graph.addEdge(normalizeGeneQuantJob, removeJob);
            graph.addEdge(normalizeISOFormQuantJob, removeJob);
        } catch (WorkflowException e1) {
            e1.printStackTrace();
        }

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
        File srcSiteResourcesImagesDir = new File("../src/site/resources/images");
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
