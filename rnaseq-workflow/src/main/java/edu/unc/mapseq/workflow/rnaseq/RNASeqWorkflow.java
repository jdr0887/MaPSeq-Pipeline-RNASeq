package edu.unc.mapseq.workflow.rnaseq;

import java.io.File;
import java.util.List;
import java.util.ResourceBundle;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.renci.jlrm.condor.CondorJob;
import org.renci.jlrm.condor.CondorJobBuilder;
import org.renci.jlrm.condor.CondorJobEdge;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.dao.model.Flowcell;
import edu.unc.mapseq.dao.model.Sample;
import edu.unc.mapseq.dao.model.WorkflowRunAttempt;
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
import edu.unc.mapseq.module.mapsplice.QualityScaleType;
import edu.unc.mapseq.module.mapsplice.RSEMCalculateExpressionCLI;
import edu.unc.mapseq.module.mapsplice.ReadChromoSizeCLI;
import edu.unc.mapseq.module.mapsplice.ReadsToUnmappedSAMCLI;
import edu.unc.mapseq.module.mapsplice.SetUnmappedBitFlagCLI;
import edu.unc.mapseq.module.picard.PicardAddOrReplaceReadGroupsCLI;
import edu.unc.mapseq.module.picard.PicardSortOrderType;
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
import edu.unc.mapseq.workflow.impl.AbstractSampleWorkflow;
import edu.unc.mapseq.workflow.impl.WorkflowJobFactory;
import edu.unc.mapseq.workflow.impl.WorkflowUtil;

public class RNASeqWorkflow extends AbstractSampleWorkflow {

    private final Logger logger = LoggerFactory.getLogger(RNASeqWorkflow.class);

    public RNASeqWorkflow() {
        super();
    }

    @Override
    public String getName() {
        return RNASeqWorkflow.class.getSimpleName().replace("Workflow", "");
    }

    @Override
    public String getVersion() {
        ResourceBundle bundle = ResourceBundle.getBundle("edu/unc/mapseq/workflow/rnaseq/workflow");
        String version = bundle.getString("version");
        return StringUtils.isNotEmpty(version) ? version : "0.0.1-SNAPSHOT";
    }

    @Override
    public void preRun() throws WorkflowException {
        super.preRun();

        Set<Sample> sampleSet = getAggregatedSamples();
        logger.info("sampleSet.size(): {}", sampleSet.size());

        for (Sample sample : sampleSet) {

            if ("Undetermined".equals(sample.getBarcode())) {
                continue;
            }

            File outputDirectory = new File(sample.getOutputDirectory(), getName());
            File tmpDirectory = new File(outputDirectory, "tmp");
            tmpDirectory.mkdirs();

            File tmpDir = new File(outputDirectory, "tmp");

            File originalDir = new File(tmpDir, "original");
            originalDir.mkdirs();

            File bestDir = new File(tmpDir, "best");
            bestDir.mkdirs();

            File remapDir = new File(tmpDir, "remap");
            remapDir.mkdirs();

            File fusionDir = new File(tmpDir, "fusion");

            File fusionDataDir = new File(fusionDir, "data");

            File fusionDataPERDir = new File(fusionDataDir, "PER");
            fusionDataPERDir.mkdirs();

            File fusionDataSingleDir = new File(fusionDataDir, "single");
            fusionDataSingleDir.mkdirs();

            File fusionResultDir = new File(fusionDir, "result");

            File fusionResultPERProbDir = new File(fusionResultDir, "PER_prob");
            fusionResultPERProbDir.mkdirs();

            File fusionResultFusionReadDir = new File(fusionResultDir, "fusionRead");
            fusionResultFusionReadDir.mkdirs();

            File fusionResultJunctionSupportDir = new File(fusionResultDir, "junction_support");
            fusionResultJunctionSupportDir.mkdirs();

            File clusterDir = new File(tmpDir, "cluster");
            clusterDir.mkdirs();

            File resultDir = new File(clusterDir, "result");

            File pairedEndMatchDir = new File(resultDir, "PEmatch");
            pairedEndMatchDir.mkdirs();

            File pairedEndMatchRegionBEDDir = new File(resultDir, "PEmatch_region_BED");
            pairedEndMatchRegionBEDDir.mkdirs();

            File pairedEndMatchRegionDir = new File(resultDir, "PEmatchregion");
            pairedEndMatchRegionDir.mkdirs();

            File pairedEndMatchRegionDetailDir = new File(resultDir, "PEmatchregion_detail");
            pairedEndMatchRegionDetailDir.mkdirs();

            File distanceDir = new File(resultDir, "distance");
            distanceDir.mkdirs();

        }

    }

    @Override
    public Graph<CondorJob, CondorJobEdge> createGraph() throws WorkflowException {
        logger.info("ENTERING createGraph()");

        DirectedGraph<CondorJob, CondorJobEdge> graph = new DefaultDirectedGraph<CondorJob, CondorJobEdge>(
                CondorJobEdge.class);

        int count = 0;

        Set<Sample> sampleSet = getAggregatedSamples();
        logger.info("sampleSet.size(): {}", sampleSet.size());

        String siteName = getWorkflowBeanService().getAttributes().get("siteName");
        String bowtieIndexDirectory = getWorkflowBeanService().getAttributes().get("bowtieIndexDirectory");
        String chromosomeDirectory = getWorkflowBeanService().getAttributes().get("chromosomeDirectory");
        String junctions = getWorkflowBeanService().getAttributes().get("junctions");
        String compositeExons = getWorkflowBeanService().getAttributes().get("compositeExons");
        String bed = getWorkflowBeanService().getAttributes().get("bed");
        String referenceSequencePrefix = getWorkflowBeanService().getAttributes().get("referenceSequencePrefix");
        String readGroupPlatform = getWorkflowBeanService().getAttributes().get("readGroupPlatform");
        String readGroupPlatformUnit = getWorkflowBeanService().getAttributes().get("readGroupPlatformUnit");

        WorkflowRunAttempt attempt = getWorkflowRunAttempt();

        for (Sample sample : sampleSet) {

            if ("Undetermined".equals(sample.getBarcode())) {
                continue;
            }

            logger.debug(sample.toString());

            Flowcell flowcell = sample.getFlowcell();
            File outputDirectory = new File(sample.getOutputDirectory(), getName());
            File tmpDirectory = new File(outputDirectory, "tmp");
            tmpDirectory.mkdirs();

            File tmpDir = new File(outputDirectory, "tmp");
            File originalDir = new File(tmpDir, "original");
            File bestDir = new File(tmpDir, "best");
            File remapDir = new File(tmpDir, "remap");
            File fusionDir = new File(tmpDir, "fusion");
            File clusterDir = new File(tmpDir, "cluster");
            File resultDir = new File(clusterDir, "result");

            List<File> readPairList = WorkflowUtil.getReadPairList(sample.getFileDatas(), flowcell.getName(),
                    sample.getLaneIndex());
            logger.info("fileList = {}", readPairList.size());

            // assumption: a dash is used as a delimiter between a participantId
            // and the external code
            int idx = sample.getName().lastIndexOf("-");
            String sampleName = idx != -1 ? sample.getName().substring(0, idx) : sample.getName();

            if (readPairList.size() != 2) {
                logger.error("readPairList.size() != 2");
                throw new WorkflowException("Need two fastq files to process...check FileData/Entity mapping");
            }

            File r1FastqFile = readPairList.get(0);
            String r1FastqRootName = WorkflowUtil.getRootFastqName(r1FastqFile.getName());

            File r2FastqFile = readPairList.get(1);
            String r2FastqRootName = WorkflowUtil.getRootFastqName(r2FastqFile.getName());

            String fastqLaneRootName = StringUtils.removeEnd(r2FastqRootName, "_R2");

            // new job
            CondorJobBuilder builder = WorkflowJobFactory.createJob(++count, GUnZipCLI.class, attempt.getId(),
                    sample.getId()).siteName(siteName);
            File gunzippedFastqR1 = new File(outputDirectory, r1FastqRootName + ".fastq");
            builder.addArgument(GUnZipCLI.GZFILE, r1FastqFile.getAbsolutePath()).addArgument(GUnZipCLI.EXTRACTFILE,
                    gunzippedFastqR1.getAbsolutePath());
            CondorJob gunzipFastqR1Job = builder.build();
            logger.info(gunzipFastqR1Job.toString());
            graph.addVertex(gunzipFastqR1Job);

            // new job
            builder = WorkflowJobFactory
                    .createJob(++count, UBUFastqFormatterCLI.class, attempt.getId(), sample.getId()).siteName(siteName);
            File fastqFormatterR1Out = new File(outputDirectory, gunzippedFastqR1.getName().replace(".fastq",
                    ".filtered.fastq"));
            builder.addArgument(UBUFastqFormatterCLI.INPUT, gunzippedFastqR1.getAbsolutePath())
                    .addArgument(UBUFastqFormatterCLI.STRIP).addArgument(UBUFastqFormatterCLI.SUFFIX, "/1")
                    .addArgument(UBUFastqFormatterCLI.OUTPUT, fastqFormatterR1Out.getAbsolutePath());
            CondorJob fastqFormatterR1Job = builder.build();
            logger.info(fastqFormatterR1Job.toString());
            graph.addVertex(fastqFormatterR1Job);
            graph.addEdge(gunzipFastqR1Job, fastqFormatterR1Job);

            // new job
            builder = WorkflowJobFactory.createJob(++count, GUnZipCLI.class, attempt.getId(), sample.getId()).siteName(
                    siteName);
            File gunzippedFastqR2 = new File(outputDirectory, r2FastqRootName + ".fastq");
            builder.addArgument(GUnZipCLI.GZFILE, r2FastqFile.getAbsolutePath()).addArgument(GUnZipCLI.EXTRACTFILE,
                    gunzippedFastqR2.getAbsolutePath());
            CondorJob gunzipFastqR2Job = builder.build();
            logger.info(gunzipFastqR2Job.toString());
            graph.addVertex(gunzipFastqR2Job);

            // new job
            builder = WorkflowJobFactory
                    .createJob(++count, UBUFastqFormatterCLI.class, attempt.getId(), sample.getId()).siteName(siteName);
            File fastqFormatterR2Out = new File(outputDirectory, gunzippedFastqR2.getName().replace(".fastq",
                    ".filtered.fastq"));
            builder.addArgument(UBUFastqFormatterCLI.INPUT, gunzippedFastqR2.getAbsolutePath())
                    .addArgument(UBUFastqFormatterCLI.STRIP).addArgument(UBUFastqFormatterCLI.SUFFIX, "/2")
                    .addArgument(UBUFastqFormatterCLI.OUTPUT, fastqFormatterR2Out.getAbsolutePath());
            CondorJob fastqFormatterR2Job = builder.build();
            logger.info(fastqFormatterR2Job.toString());
            graph.addVertex(fastqFormatterR2Job);
            graph.addEdge(gunzipFastqR2Job, fastqFormatterR2Job);

            // new job
            builder = WorkflowJobFactory.createJob(++count, DetermineReadLengthCLI.class, attempt.getId(),
                    sample.getId()).siteName(siteName);
            File determineReadLengthOutput = new File(tmpDir, "readLengthProps.xml");
            builder.addArgument(DetermineReadLengthCLI.INPUT, fastqFormatterR1Out.getAbsolutePath())
                    .addArgument(DetermineReadLengthCLI.INPUT, fastqFormatterR2Out.getAbsolutePath())
                    .addArgument(DetermineReadLengthCLI.OUTPUT, determineReadLengthOutput.getAbsolutePath());
            CondorJob determineReadLengthJob = builder.build();
            logger.info(determineReadLengthJob.toString());
            graph.addVertex(determineReadLengthJob);
            graph.addEdge(fastqFormatterR1Job, determineReadLengthJob);
            graph.addEdge(fastqFormatterR2Job, determineReadLengthJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, ReadChromoSizeCLI.class, attempt.getId(), sample.getId())
                    .siteName(siteName);
            File chromosomeHeadOutput = new File(tmpDir, "chrom_head");
            File chromosomeIndexOutput = new File(tmpDir, "chromo.fai");
            File chromosomeNamesOutput = new File(fusionDir, "chrName.txt");
            File chromosomeSizesOutput = new File(tmpDir, "chrom_sizes");
            builder.addArgument(ReadChromoSizeCLI.CHROMOSOMEHEADOUTPUT, chromosomeHeadOutput.getAbsolutePath())
                    .addArgument(ReadChromoSizeCLI.CHROMOSOMEINDEX, chromosomeIndexOutput.getAbsolutePath())
                    .addArgument(ReadChromoSizeCLI.CHROMOSOMENAMESOUTPUT, chromosomeNamesOutput.getAbsolutePath())
                    .addArgument(ReadChromoSizeCLI.CHROMOSOMESIZESOUTPUT, chromosomeSizesOutput.getAbsolutePath());
            File chromDir = new File(chromosomeDirectory);
            for (File f : chromDir.listFiles()) {
                if (f.getName().endsWith(".fa")) {
                    builder.addArgument(ReadChromoSizeCLI.SAMFILE, f.getAbsolutePath());
                }
            }
            CondorJob readChromoSizesJob = builder.build();
            logger.info(readChromoSizesJob.toString());
            graph.addVertex(readChromoSizesJob);

            // new job
            builder = WorkflowJobFactory
                    .createJob(++count, MapSpliceMultiThreadCLI.class, attempt.getId(), sample.getId())
                    .siteName(siteName).numberOfProcessors(8);
            File mapspliceMultiThreadJobOutput = new File(originalDir, "bowtie.output");
            File originalSAMMapspliceMultiThreadOutput = new File(originalDir, "original.sam");
            builder.addArgument(MapSpliceMultiThreadCLI.FASTQ1, fastqFormatterR1Out.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.FASTQ2, fastqFormatterR2Out.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.THREADS, 8)
                    .addArgument(MapSpliceMultiThreadCLI.QUALITYSCALE, QualityScaleType.phred33)
                    .addArgument(MapSpliceMultiThreadCLI.MINIMUMINTRON, 50)
                    .addArgument(MapSpliceMultiThreadCLI.SUPPRESSALIGNMENTSOVER, 40)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMINTRONSINGLE, 200000)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMINTRONDOUBLE, 200000)
                    .addArgument(MapSpliceMultiThreadCLI.MINIMUMLENGTH, 25)
                    .addArgument(MapSpliceMultiThreadCLI.SPLICEMISMATCHES, 1)
                    .addArgument(MapSpliceMultiThreadCLI.SPLICEONLY)
                    .addArgument(MapSpliceMultiThreadCLI.REPORTUPTONALIGNMENTSPERREAD, 40)
                    .addArgument(MapSpliceMultiThreadCLI.FASTQREADFORMAT)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMINSERTIONS, 6)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMDELETIONS, 6)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMAPPENDMISMATCHES, 3)
                    .addArgument(MapSpliceMultiThreadCLI.SEGMENTLENGTH, 25)
                    .addArgument(MapSpliceMultiThreadCLI.OUTPUT, mapspliceMultiThreadJobOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.MAPSPLICEOUT,
                            originalSAMMapspliceMultiThreadOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.REFERENCESEQUENCEDIRECTORY, chromosomeDirectory)
                    .addArgument(MapSpliceMultiThreadCLI.CHROMOSOMESIZE, chromosomeSizesOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.INDEX, bowtieIndexDirectory);
            CondorJob originalSAMMapspliceMultiThreadJob = builder.build();
            logger.info(originalSAMMapspliceMultiThreadJob.toString());
            graph.addVertex(originalSAMMapspliceMultiThreadJob);
            graph.addEdge(fastqFormatterR1Job, originalSAMMapspliceMultiThreadJob);
            graph.addEdge(fastqFormatterR2Job, originalSAMMapspliceMultiThreadJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, RegexCatCLI.class, attempt.getId()).siteName(siteName);
            builder.addArgument(RegexCatCLI.DIRECTORY, originalDir.getAbsolutePath())
                    .addArgument(RegexCatCLI.OUTPUT, originalSAMMapspliceMultiThreadOutput.getAbsolutePath())
                    .addArgument(RegexCatCLI.REGEX, "^.*\\.sam\\.[0-9]");
            CondorJob originalSAMRegexCatJob = builder.build();
            logger.info(originalSAMRegexCatJob.toString());
            graph.addVertex(originalSAMRegexCatJob);
            graph.addEdge(originalSAMMapspliceMultiThreadJob, originalSAMRegexCatJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, NewSAM2JuntionCLI.class, attempt.getId())
                    .siteName(siteName);
            File originalAllJunctionsFile = new File(originalDir, "ori.all_junctions.txt");
            builder.addArgument(NewSAM2JuntionCLI.CHROMOSOMEFILESDIRECTORY, chromosomeDirectory)
                    .addArgument(NewSAM2JuntionCLI.INPUT, originalSAMMapspliceMultiThreadOutput.getAbsolutePath())
                    .addArgument(NewSAM2JuntionCLI.JUNCTIONFILE, originalAllJunctionsFile.getAbsolutePath())
                    .addArgument(NewSAM2JuntionCLI.MINIMUMANCHOR, 1).addArgument(NewSAM2JuntionCLI.MINIMUMINTRON, 1)
                    .addArgument(NewSAM2JuntionCLI.MAXIMUMINTRON, 350000000)
                    .addArgument(NewSAM2JuntionCLI.READLENGTH, determineReadLengthOutput.getAbsolutePath());
            CondorJob sam2JunctionArrayJob = builder.build();
            logger.info(sam2JunctionArrayJob.toString());
            graph.addVertex(sam2JunctionArrayJob);
            graph.addEdge(determineReadLengthJob, sam2JunctionArrayJob);
            graph.addEdge(originalSAMRegexCatJob, sam2JunctionArrayJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, Filter1HitsCLI.class, attempt.getId()).siteName(siteName);
            File canonicalFilter1HitsOutput = new File(bestDir,
                    "ori.all_junctions.filtered_by_min_mis_lpq.remained.txt");
            File nonCanonicalFilter1HitsOutput = new File(bestDir,
                    "ori.all_junctions.filtered_by_min_mis_lpq.filtered.txt");
            builder.addArgument(Filter1HitsCLI.JUNCTIONFILE, originalAllJunctionsFile.getAbsolutePath())
                    .addArgument(Filter1HitsCLI.CANONICALJUNCTIONOUTPUT, canonicalFilter1HitsOutput.getAbsolutePath())
                    .addArgument(Filter1HitsCLI.NONCANONICALJUNCTIONOUTPUT,
                            nonCanonicalFilter1HitsOutput.getAbsolutePath())
                    .addArgument(Filter1HitsCLI.MINIMUMLPQ, 0.3).addArgument(Filter1HitsCLI.MINIMUMMISMATCH, 2);
            CondorJob filter1HitsJob = builder.build();
            logger.info(filter1HitsJob.toString());
            graph.addVertex(filter1HitsJob);
            graph.addEdge(sam2JunctionArrayJob, filter1HitsJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, FilterJunctionByROCarguNonCanonicalCLI.class,
                    attempt.getId()).siteName(siteName);
            File bestJunctionOutput = new File(bestDir, "best_junction.txt");
            File nonCanonicalBestJunctionOutput = new File(bestDir,
                    "best_junction_semi_non_canon_filtered_by_ROCargu.txt");
            builder.addArgument(FilterJunctionByROCarguNonCanonicalCLI.ENTROPYWEIGHT, 0.097718)
                    .addArgument(FilterJunctionByROCarguNonCanonicalCLI.LPQWEIGHT, 0.66478)
                    .addArgument(FilterJunctionByROCarguNonCanonicalCLI.MINSCORE, 0.719)
                    .addArgument(FilterJunctionByROCarguNonCanonicalCLI.AVEMISMATCHWEIGHT, -0.21077)
                    .addArgument(FilterJunctionByROCarguNonCanonicalCLI.JUNCTIONFILE,
                            canonicalFilter1HitsOutput.getAbsolutePath())
                    .addArgument(FilterJunctionByROCarguNonCanonicalCLI.CANONICALJUNCTIONOUTPUT,
                            bestJunctionOutput.getAbsolutePath())
                    .addArgument(FilterJunctionByROCarguNonCanonicalCLI.NONCANONICALJUNCTIONOUTPUT,
                            nonCanonicalBestJunctionOutput.getAbsolutePath())
                    .addArgument(FilterJunctionByROCarguNonCanonicalCLI.INTRONWEIGHT, 0)
                    .addArgument(FilterJunctionByROCarguNonCanonicalCLI.SUMLENGTHWEIGHT, 0)
                    .addArgument(FilterJunctionByROCarguNonCanonicalCLI.MINFLANKCASE, 5);
            CondorJob filterJunctionByROCarguNonCanonicalJob = builder.build();
            logger.info(filterJunctionByROCarguNonCanonicalJob.toString());
            graph.addVertex(filterJunctionByROCarguNonCanonicalJob);
            graph.addEdge(filter1HitsJob, filterJunctionByROCarguNonCanonicalJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, JunctionSequenceConstructionCLI.class, attempt.getId())
                    .siteName(siteName);
            File junctionSequenceConstructionOutput = new File(remapDir, "synthetic_alljunc_sequence.txt");
            builder.addArgument(JunctionSequenceConstructionCLI.JUNCTION, bestJunctionOutput.getAbsolutePath())
                    .addArgument(JunctionSequenceConstructionCLI.REFERENCESEQUENCEDIRECTORY, chromosomeDirectory)
                    .addArgument(JunctionSequenceConstructionCLI.OUTPUT,
                            junctionSequenceConstructionOutput.getAbsolutePath())
                    .addArgument(JunctionSequenceConstructionCLI.MINIMUMANCHOR, 2)
                    .addArgument(JunctionSequenceConstructionCLI.MAXIMUMANCHOR, 38)
                    .addArgument(JunctionSequenceConstructionCLI.MAXIMUMSEQUENCETHRESHOLD, 400);
            CondorJob junctionSequenceConstructionJob = builder.build();
            logger.info(junctionSequenceConstructionJob.toString());
            graph.addVertex(junctionSequenceConstructionJob);
            graph.addEdge(filterJunctionByROCarguNonCanonicalJob, junctionSequenceConstructionJob);

            builder = WorkflowJobFactory.createJob(++count, BowtieBuildCLI.class, attempt.getId()).siteName(siteName);
            File synIndexPrefix = new File(remapDir, "syn_idx_prefix");
            builder.addArgument(BowtieBuildCLI.INPUT, junctionSequenceConstructionOutput.getAbsolutePath())
                    .addArgument(BowtieBuildCLI.PREFIX, synIndexPrefix.getAbsolutePath());
            CondorJob bowtieBuildJob = builder.build();
            logger.info(bowtieBuildJob.toString());
            graph.addVertex(bowtieBuildJob);
            graph.addEdge(junctionSequenceConstructionJob, bowtieBuildJob);

            // new job
            builder = WorkflowJobFactory
                    .createJob(++count, MapSpliceMultiThreadCLI.class, attempt.getId(), sample.getId())
                    .siteName(siteName).numberOfProcessors(8);
            File unmappedOutput = new File(remapDir, "remap_unmapped");
            mapspliceMultiThreadJobOutput = new File(remapDir, "bowtie.output");
            File remappedSAMMapspliceMultiThreadOutput = new File(remapDir, "remapped.sam");
            builder.addArgument(MapSpliceMultiThreadCLI.FASTQ1, fastqFormatterR1Out.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.FASTQ2, fastqFormatterR2Out.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.THREADS, 8)
                    .addArgument(MapSpliceMultiThreadCLI.QUALITYSCALE, QualityScaleType.phred33)
                    .addArgument(MapSpliceMultiThreadCLI.OUTPUTUNMAPPED, unmappedOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.MINIMUMINTRON, 50)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMINTRONSINGLE, 200000)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMINTRONDOUBLE, 200000)
                    .addArgument(MapSpliceMultiThreadCLI.MINIMUMLENGTH, 25)
                    .addArgument(MapSpliceMultiThreadCLI.SEGMENTLENGTH, 25)
                    .addArgument(MapSpliceMultiThreadCLI.SUPPRESSALIGNMENTSOVER, 40)
                    .addArgument(MapSpliceMultiThreadCLI.REPORTUPTONALIGNMENTSPERREAD, 40)
                    .addArgument(MapSpliceMultiThreadCLI.OPTIMIZEREPEATS)
                    .addArgument(MapSpliceMultiThreadCLI.JUNCTIONINDEX, synIndexPrefix.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.SPLICEMISMATCHES, 1)
                    .addArgument(MapSpliceMultiThreadCLI.FASTQREADFORMAT)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMINSERTIONS, 6)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMDELETIONS, 6)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMAPPENDMISMATCHES, 3)
                    .addArgument(MapSpliceMultiThreadCLI.OUTPUT, mapspliceMultiThreadJobOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.MAPSPLICEOUT,
                            remappedSAMMapspliceMultiThreadOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.REFERENCESEQUENCEDIRECTORY, chromosomeDirectory)
                    .addArgument(MapSpliceMultiThreadCLI.CHROMOSOMESIZE, chromosomeSizesOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.INDEX, bowtieIndexDirectory);
            CondorJob remappedSAMMapspliceMultiThreadJob = builder.build();
            logger.info(remappedSAMMapspliceMultiThreadJob.toString());
            graph.addVertex(remappedSAMMapspliceMultiThreadJob);
            graph.addEdge(bowtieBuildJob, remappedSAMMapspliceMultiThreadJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, RegexCatCLI.class, attempt.getId()).siteName(siteName);
            builder.addArgument(RegexCatCLI.DIRECTORY, remapDir.getAbsolutePath())
                    .addArgument(RegexCatCLI.OUTPUT, remappedSAMMapspliceMultiThreadOutput.getAbsolutePath())
                    .addArgument(RegexCatCLI.REGEX, "^remapped\\.sam\\.[0-9]");
            CondorJob remappedRegexCatJob = builder.build();
            logger.info(remappedRegexCatJob.toString());
            graph.addVertex(remappedRegexCatJob);
            graph.addEdge(remappedSAMMapspliceMultiThreadJob, remappedRegexCatJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, RegexCatCLI.class, attempt.getId()).siteName(siteName);
            File remapUnmapped1RegexCatOutput = new File(remapDir, "remap_unmapped.1");
            builder.addArgument(RegexCatCLI.DIRECTORY, remapDir.getAbsolutePath())
                    .addArgument(RegexCatCLI.OUTPUT, remapUnmapped1RegexCatOutput.getAbsolutePath())
                    .addArgument(RegexCatCLI.REGEX, "^remap_unmapped\\.[0-9]\\.1");
            CondorJob remapUnmapped1RegexCatJob = builder.build();
            logger.info(remapUnmapped1RegexCatJob.toString());
            graph.addVertex(remapUnmapped1RegexCatJob);
            graph.addEdge(remappedSAMMapspliceMultiThreadJob, remapUnmapped1RegexCatJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, RegexCatCLI.class, attempt.getId()).siteName(siteName);
            File remapUnmapped2RegexCatOutput = new File(remapDir, "remap_unmapped.2");
            builder.addArgument(RegexCatCLI.DIRECTORY, remapDir.getAbsolutePath())
                    .addArgument(RegexCatCLI.OUTPUT, remapUnmapped2RegexCatOutput.getAbsolutePath())
                    .addArgument(RegexCatCLI.REGEX, "^remap_unmapped\\.[0-9]\\.2");
            CondorJob remapUnmapped2RegexCatJob = builder.build();
            logger.info(remapUnmapped2RegexCatJob.toString());
            graph.addVertex(remapUnmapped2RegexCatJob);
            graph.addEdge(remappedSAMMapspliceMultiThreadJob, remapUnmapped2RegexCatJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, AlignmentHandlerMultiCLI.class, attempt.getId())
                    .siteName(siteName).numberOfProcessors(8);
            File filteredAlignmentBase = new File(remapDir, "_filtered_normal_alignments");
            builder.addArgument(AlignmentHandlerMultiCLI.ADDSOFTCLIP, 1)
                    .addArgument(AlignmentHandlerMultiCLI.AVERAGEFRAGMENTLENGTH, 225)
                    .addArgument(AlignmentHandlerMultiCLI.BOUNDARY, 36)
                    .addArgument(AlignmentHandlerMultiCLI.CHROMOSOMESIZEFILE, chromosomeSizesOutput.getAbsolutePath())
                    .addArgument(AlignmentHandlerMultiCLI.ENCOMPASSINGFUSIONREGIONEXTENSION, 50000)
                    .addArgument(AlignmentHandlerMultiCLI.FILTEREDALIGNMENTAPPEND, "_filtered_normal_alignments")
                    .addArgument(AlignmentHandlerMultiCLI.FILTEREDALIGNMENTBASE,
                            filteredAlignmentBase.getAbsolutePath())
                    .addArgument(AlignmentHandlerMultiCLI.FILTERFLAG, 12 + 32 + 256)
                    .addArgument(AlignmentHandlerMultiCLI.FRAGMENTLENGTH, 400)
                    .addArgument(AlignmentHandlerMultiCLI.FRAGMENTLENGTHSD, 100)
                    .addArgument(AlignmentHandlerMultiCLI.INPUT,
                            remappedSAMMapspliceMultiThreadOutput.getAbsolutePath())
                    .addArgument(AlignmentHandlerMultiCLI.INTRONDISTANCESD, 500)
                    .addArgument(AlignmentHandlerMultiCLI.MATEDISTANCESD, 100)
                    .addArgument(AlignmentHandlerMultiCLI.MAXIMUMANCHORDIFF, 50)
                    .addArgument(AlignmentHandlerMultiCLI.MAXIMUMDELETION, 6)
                    .addArgument(AlignmentHandlerMultiCLI.MAXIMUMHITS, 40)
                    .addArgument(AlignmentHandlerMultiCLI.MAXIMUMMATEDISTANCE, 50000)
                    .addArgument(AlignmentHandlerMultiCLI.READLENGTHPROPERTIES,
                            determineReadLengthOutput.getAbsolutePath())
                    .addArgument(AlignmentHandlerMultiCLI.MINIMUMANCHOR, 0)
                    .addArgument(AlignmentHandlerMultiCLI.MINIMUMCOVERAGE, 0)
                    .addArgument(AlignmentHandlerMultiCLI.MINIMUMENCOMPASSCOUNT, 1)
                    .addArgument(AlignmentHandlerMultiCLI.MINIMUMENTROPY, -0.0001)
                    .addArgument(AlignmentHandlerMultiCLI.MINIMUMINSERTION, 6)
                    .addArgument(AlignmentHandlerMultiCLI.MINIMUMJUNCTIONANCHOR, 10)
                    .addArgument(AlignmentHandlerMultiCLI.MINIMUMMISMATCH, 5)
                    .addArgument(AlignmentHandlerMultiCLI.THREADS, 8);
            CondorJob alignmentHandlerMultiJob = builder.build();
            logger.info(alignmentHandlerMultiJob.toString());
            graph.addVertex(alignmentHandlerMultiJob);
            graph.addEdge(determineReadLengthJob, alignmentHandlerMultiJob);
            graph.addEdge(remappedRegexCatJob, alignmentHandlerMultiJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, CatCLI.class, attempt.getId()).siteName(siteName);
            File fusionPaired = new File(remapDir, remappedSAMMapspliceMultiThreadOutput.getName()
                    + "_filtered_normal_alignments.fusion_paired");
            File bothUnsplicedFusionPairedAlignments = new File(remapDir,
                    remappedSAMMapspliceMultiThreadOutput.getName()
                            + "_filtered_normal_alignments.bothunspliced.fusion_paired");
            File fusionPairedAlignments = new File(remapDir, remappedSAMMapspliceMultiThreadOutput.getName()
                    + "_filtered_normal_alignments.fusion_paired.comb");
            builder.addArgument(CatCLI.FILES, fusionPaired.getAbsolutePath())
                    .addArgument(CatCLI.FILES, bothUnsplicedFusionPairedAlignments.getAbsolutePath())
                    .addArgument(CatCLI.OUTPUT, fusionPairedAlignments.getAbsolutePath());
            CondorJob filteredNormalAlignmentsFusionPairedCatJob = builder.build();
            logger.info(filteredNormalAlignmentsFusionPairedCatJob.toString());
            graph.addVertex(filteredNormalAlignmentsFusionPairedCatJob);
            graph.addEdge(alignmentHandlerMultiJob, filteredNormalAlignmentsFusionPairedCatJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, CatCLI.class, attempt.getId()).siteName(siteName);
            File single = new File(remapDir, remappedSAMMapspliceMultiThreadOutput.getName()
                    + "_filtered_normal_alignments.single");
            File bothUnsplicedSingleAlignments = new File(remapDir, remappedSAMMapspliceMultiThreadOutput.getName()
                    + "_filtered_normal_alignments.bothunspliced.single");
            File singleAlignments = new File(remapDir, remappedSAMMapspliceMultiThreadOutput.getName()
                    + "_filtered_normal_alignments.single.comb");
            builder.addArgument(CatCLI.FILES, single.getAbsolutePath())
                    .addArgument(CatCLI.FILES, bothUnsplicedSingleAlignments.getAbsolutePath())
                    .addArgument(CatCLI.OUTPUT, singleAlignments.getAbsolutePath());
            CondorJob filteredNormalAlignmentsSingleCatJob = builder.build();
            logger.info(filteredNormalAlignmentsSingleCatJob.toString());
            graph.addVertex(filteredNormalAlignmentsSingleCatJob);
            graph.addEdge(alignmentHandlerMultiJob, filteredNormalAlignmentsSingleCatJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, CatCLI.class, attempt.getId()).siteName(siteName);
            File paired = new File(remapDir, remappedSAMMapspliceMultiThreadOutput.getName()
                    + "_filtered_normal_alignments.paired");
            File bothUnsplicedPairedAlignments = new File(remapDir, remappedSAMMapspliceMultiThreadOutput.getName()
                    + "_filtered_normal_alignments.bothunspliced.paired");
            File pairedAlignments = new File(remapDir, remappedSAMMapspliceMultiThreadOutput.getName()
                    + "_filtered_normal_alignments.paired.comb");
            builder.addArgument(CatCLI.FILES, paired.getAbsolutePath())
                    .addArgument(CatCLI.FILES, bothUnsplicedPairedAlignments.getAbsolutePath())
                    .addArgument(CatCLI.OUTPUT, pairedAlignments.getAbsolutePath());
            CondorJob filteredNormalAlignmentsPairedCatJob = builder.build();
            logger.info(filteredNormalAlignmentsPairedCatJob.toString());
            graph.addVertex(filteredNormalAlignmentsPairedCatJob);
            graph.addEdge(alignmentHandlerMultiJob, filteredNormalAlignmentsPairedCatJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, SedCLI.class, attempt.getId()).siteName(siteName);
            File fusionSAM = new File(fusionDir, "sed.fusion.sam");
            builder.addArgument(SedCLI.SOURCE, fusionPairedAlignments.getAbsolutePath())
                    .addArgument(SedCLI.REGULAREXPRESSION, "s/^/1~/")
                    .addArgument(SedCLI.OUTPUT, fusionSAM.getAbsolutePath());
            CondorJob sedJob = builder.build();
            logger.info(sedJob.toString());
            graph.addVertex(sedJob);
            graph.addEdge(filteredNormalAlignmentsFusionPairedCatJob, sedJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, ParseClusterCLI.class, attempt.getId()).siteName(siteName);
            builder.addArgument(ParseClusterCLI.INPUT, fusionSAM.getAbsolutePath()).addArgument(
                    ParseClusterCLI.OUTPUTDIRECTORY, clusterDir.getAbsolutePath());
            CondorJob parseClusterJob = builder.build();
            logger.info(parseClusterJob.toString());
            graph.addVertex(parseClusterJob);
            graph.addEdge(sedJob, parseClusterJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, ClusterCLI.class, attempt.getId()).siteName(siteName);
            builder.addArgument(ClusterCLI.CLUSTERDIRECTORY, clusterDir.getAbsolutePath());
            CondorJob clusterJob = builder.build();
            logger.info(clusterJob.toString());
            graph.addVertex(clusterJob);
            graph.addEdge(parseClusterJob, clusterJob);

            // new job
            builder = WorkflowJobFactory
                    .createJob(++count, MapSpliceMultiThreadCLI.class, attempt.getId(), sample.getId())
                    .siteName(siteName).numberOfProcessors(8);
            File normalSAMMapspliceMultiThreadFusionOutput = new File(fusionDir, "normal.sam");
            unmappedOutput = new File(fusionDir, "fusion_original_unmapped");
            mapspliceMultiThreadJobOutput = new File(fusionDir, "bowtie_original.output");
            File clusterRegion = new File(resultDir, "cluster.txt");
            File fusion = new File(fusionDir, "fusion_alignments_original.sam");
            builder.addArgument(MapSpliceMultiThreadCLI.FASTQ1, remapUnmapped1RegexCatOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.FASTQ2, remapUnmapped2RegexCatOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.THREADS, 8)
                    .addArgument(MapSpliceMultiThreadCLI.QUALITYSCALE, QualityScaleType.phred33)
                    .addArgument(MapSpliceMultiThreadCLI.MAPSPLICEOUT,
                            normalSAMMapspliceMultiThreadFusionOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.OUTPUTUNMAPPED, unmappedOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.OUTPUT, mapspliceMultiThreadJobOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.OPTIMIZEREPEATS)
                    .addArgument(MapSpliceMultiThreadCLI.CLUSTERREGION, clusterRegion.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.FUSION, fusion.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.JUNCTIONINDEX, synIndexPrefix.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.MINIMUMINTRON, 50)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMINTRONSINGLE, 200000)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMINTRONDOUBLE, 200000)
                    .addArgument(MapSpliceMultiThreadCLI.MINIMUMLENGTH, 25)
                    .addArgument(MapSpliceMultiThreadCLI.SEGMENTLENGTH, 25)
                    .addArgument(MapSpliceMultiThreadCLI.FASTQREADFORMAT)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMINSERTIONS, 6)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMDELETIONS, 6)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMAPPENDMISMATCHES, 3)
                    .addArgument(MapSpliceMultiThreadCLI.REFERENCESEQUENCEDIRECTORY, chromosomeDirectory)
                    .addArgument(MapSpliceMultiThreadCLI.REPORTUPTONALIGNMENTSPERREAD, 40)
                    .addArgument(MapSpliceMultiThreadCLI.SUPPRESSALIGNMENTSOVER, 40)
                    .addArgument(MapSpliceMultiThreadCLI.SPLICEMISMATCHES, 1)
                    .addArgument(MapSpliceMultiThreadCLI.CHROMOSOMESIZE, chromosomeSizesOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.INDEX, bowtieIndexDirectory);
            CondorJob mapspliceMultiThreadFusionJob = builder.build();
            logger.info(mapspliceMultiThreadFusionJob.toString());
            graph.addVertex(mapspliceMultiThreadFusionJob);
            graph.addEdge(clusterJob, mapspliceMultiThreadFusionJob);
            graph.addEdge(remapUnmapped1RegexCatJob, mapspliceMultiThreadFusionJob);
            graph.addEdge(remapUnmapped2RegexCatJob, mapspliceMultiThreadFusionJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, FusionSAM2JunctionFilterAnchorNewFormatCLI.class,
                    attempt.getId(), sample.getId()).siteName(siteName);
            File fusionSAM2JunctionFilterAnchorNewFormatOutput = new File(fusionDir, "original_fusion_junction.txt");
            builder.addArgument(FusionSAM2JunctionFilterAnchorNewFormatCLI.JUNCTION,
                    fusionSAM2JunctionFilterAnchorNewFormatOutput.getAbsolutePath())
                    .addArgument(FusionSAM2JunctionFilterAnchorNewFormatCLI.SAM, fusion.getAbsolutePath())
                    .addArgument(FusionSAM2JunctionFilterAnchorNewFormatCLI.MINIMUMANCHOR, 1)
                    .addArgument(FusionSAM2JunctionFilterAnchorNewFormatCLI.READLENGTH,
                            determineReadLengthOutput.getAbsolutePath())
                    .addArgument(FusionSAM2JunctionFilterAnchorNewFormatCLI.REFERENCESEQUENCEDIRECTORY,
                            chromosomeDirectory);
            CondorJob fusionSAM2JunctionFilterAnchorNewFormatJob = builder.build();
            logger.info(fusionSAM2JunctionFilterAnchorNewFormatJob.toString());
            graph.addVertex(fusionSAM2JunctionFilterAnchorNewFormatJob);
            graph.addEdge(mapspliceMultiThreadFusionJob, fusionSAM2JunctionFilterAnchorNewFormatJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, FilterOriginalFusionCLI.class, attempt.getId(),
                    sample.getId()).siteName(siteName);
            File filteredJunctions = new File(fusionDir, "original_fusion_junction.filtered.txt");
            File remainingJunctions = new File(fusionDir, "original_fusion_junction.remained.txt");
            builder.addArgument(FilterOriginalFusionCLI.JUNCTION,
                    fusionSAM2JunctionFilterAnchorNewFormatOutput.getAbsolutePath())
                    .addArgument(FilterOriginalFusionCLI.MINIMUMLPQ, 0.3)
                    .addArgument(FilterOriginalFusionCLI.MINIMUMMISMATCHES, 2)
                    .addArgument(FilterOriginalFusionCLI.FILTEREDJUNCTIONS, filteredJunctions.getAbsolutePath())
                    .addArgument(FilterOriginalFusionCLI.REMAININGJUNCTIONS, remainingJunctions.getAbsolutePath());
            CondorJob filterOriginalFusionJob = builder.build();
            logger.info(filterOriginalFusionJob.toString());
            graph.addVertex(filterOriginalFusionJob);
            graph.addEdge(fusionSAM2JunctionFilterAnchorNewFormatJob, filterOriginalFusionJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, JunctionDBFusionCLI.class, attempt.getId(), sample.getId())
                    .siteName(siteName);
            File junctionDBFusionOutput = new File(fusionDir, "fusion_synthetic_sequence.txt");
            builder.addArgument(JunctionDBFusionCLI.FUSIONJUNCTION, remainingJunctions.getAbsolutePath())
                    .addArgument(JunctionDBFusionCLI.JUNCTION, bestJunctionOutput.getAbsolutePath())
                    .addArgument(JunctionDBFusionCLI.MAXIMUMANCHOR, 38)
                    .addArgument(JunctionDBFusionCLI.MAXIMUMTHRESHOLDEACH, 20)
                    .addArgument(JunctionDBFusionCLI.MAXIMUMTHRESHOLDTOTAL, 50)
                    .addArgument(JunctionDBFusionCLI.MINIMUMANCHORWIDTH, 2)
                    .addArgument(JunctionDBFusionCLI.OUTPUT, junctionDBFusionOutput.getAbsolutePath())
                    .addArgument(JunctionDBFusionCLI.REFERENCESEQUENCEDIRECTORY, chromosomeDirectory);
            CondorJob junctionDBFusionJob = builder.build();
            logger.info(junctionDBFusionJob.toString());
            graph.addVertex(junctionDBFusionJob);
            graph.addEdge(filterOriginalFusionJob, junctionDBFusionJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, BowtieBuildCLI.class, attempt.getId()).siteName(siteName);
            File synFusionIndexPrefix = new File(fusionDir, "syn_fusion_idx_prefix");
            builder.addArgument(BowtieBuildCLI.INPUT, junctionDBFusionOutput.getAbsolutePath()).addArgument(
                    BowtieBuildCLI.PREFIX, synFusionIndexPrefix.getAbsolutePath());
            CondorJob synFusionIndexBowtieBuildJob = builder.build();
            logger.info(synFusionIndexBowtieBuildJob.toString());
            graph.addVertex(synFusionIndexBowtieBuildJob);
            graph.addEdge(junctionDBFusionJob, synFusionIndexBowtieBuildJob);

            // new job
            builder = WorkflowJobFactory
                    .createJob(++count, MapSpliceMultiThreadCLI.class, attempt.getId(), sample.getId())
                    .siteName(siteName).numberOfProcessors(8);
            normalSAMMapspliceMultiThreadFusionOutput = new File(fusionDir, "normal.sam");
            unmappedOutput = new File(fusionDir, "fusion_unmapped");
            mapspliceMultiThreadJobOutput = new File(fusionDir, "bowtie.output");
            fusion = new File(fusionDir, "fusion_alignments.sam");
            builder.addArgument(MapSpliceMultiThreadCLI.FASTQ1, remapUnmapped1RegexCatOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.FASTQ2, remapUnmapped2RegexCatOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.THREADS, 8)
                    .addArgument(MapSpliceMultiThreadCLI.QUALITYSCALE, QualityScaleType.phred33)
                    .addArgument(MapSpliceMultiThreadCLI.MAPSPLICEOUT,
                            normalSAMMapspliceMultiThreadFusionOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.OUTPUTUNMAPPED, unmappedOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.OUTPUT, mapspliceMultiThreadJobOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.OPTIMIZEREPEATS)
                    .addArgument(MapSpliceMultiThreadCLI.MINIMUMINTRON, 50)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMINTRONSINGLE, 200000)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMINTRONDOUBLE, 200000)
                    .addArgument(MapSpliceMultiThreadCLI.MINIMUMLENGTH, 25)
                    .addArgument(MapSpliceMultiThreadCLI.SEGMENTLENGTH, 25)
                    .addArgument(MapSpliceMultiThreadCLI.FASTQREADFORMAT)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMINSERTIONS, 6)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMDELETIONS, 6)
                    .addArgument(MapSpliceMultiThreadCLI.MAXIMUMAPPENDMISMATCHES, 3)
                    .addArgument(MapSpliceMultiThreadCLI.REFERENCESEQUENCEDIRECTORY, chromosomeDirectory)
                    .addArgument(MapSpliceMultiThreadCLI.REPORTUPTONALIGNMENTSPERREAD, 40)
                    .addArgument(MapSpliceMultiThreadCLI.SUPPRESSALIGNMENTSOVER, 40)
                    .addArgument(MapSpliceMultiThreadCLI.SPLICEMISMATCHES, 1)
                    .addArgument(MapSpliceMultiThreadCLI.CHROMOSOMESIZE, chromosomeSizesOutput.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.CLUSTERREGION, clusterRegion.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.JUNCTIONINDEX, synIndexPrefix.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.FUSIONINDEX, synFusionIndexPrefix.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.FUSION, fusion.getAbsolutePath())
                    .addArgument(MapSpliceMultiThreadCLI.INDEX, bowtieIndexDirectory);
            mapspliceMultiThreadFusionJob = builder.build();
            logger.info(mapspliceMultiThreadFusionJob.toString());
            graph.addVertex(mapspliceMultiThreadFusionJob);
            graph.addEdge(synFusionIndexBowtieBuildJob, mapspliceMultiThreadFusionJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, SortCLI.class, attempt.getId()).siteName(siteName);
            File sortOutput = new File(fusionDir, "combined_fusion_normal.sam");
            builder.addArgument(SortCLI.OUTPUT, sortOutput.getAbsolutePath()).addArgument(SortCLI.BUFFERSIZE, 3500000)
                    .addArgument(SortCLI.KEY, "1,1").addArgument(SortCLI.TMPDIRECTORY, tmpDir.getAbsolutePath())
                    .addArgument(SortCLI.INPUT, singleAlignments.getAbsolutePath())
                    .addArgument(SortCLI.INPUT, fusionPairedAlignments.getAbsolutePath())
                    .addArgument(SortCLI.INPUT, fusion.getAbsolutePath());
            CondorJob sortJob = builder.build();
            logger.info(sortJob.toString());
            graph.addVertex(sortJob);
            graph.addEdge(mapspliceMultiThreadFusionJob, sortJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, RegexCatCLI.class, attempt.getId()).siteName(siteName);
            File fusionUnmapped1RegexCatOutput = new File(fusionDir, "fusion_unmapped.1");
            builder.addArgument(RegexCatCLI.DIRECTORY, fusionDir.getAbsolutePath())
                    .addArgument(RegexCatCLI.OUTPUT, fusionUnmapped1RegexCatOutput.getAbsolutePath())
                    .addArgument(RegexCatCLI.REGEX, "^fusion_unmapped\\.[0-9]\\.1");
            CondorJob fusionUnmapped1RegexCatJob = builder.build();
            logger.info(fusionUnmapped1RegexCatJob.toString());
            graph.addVertex(fusionUnmapped1RegexCatJob);
            graph.addEdge(mapspliceMultiThreadFusionJob, fusionUnmapped1RegexCatJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, RegexCatCLI.class, attempt.getId()).siteName(siteName);
            File fusionUnmapped2RegexCatOutput = new File(fusionDir, "fusion_unmapped.2");
            builder.addArgument(RegexCatCLI.DIRECTORY, fusionDir.getAbsolutePath())
                    .addArgument(RegexCatCLI.OUTPUT, fusionUnmapped2RegexCatOutput.getAbsolutePath())
                    .addArgument(RegexCatCLI.REGEX, "^fusion_unmapped\\.[0-9]\\.2");
            CondorJob fusionUnmapped2RegexCatJob = builder.build();
            logger.info(fusionUnmapped2RegexCatJob.toString());
            graph.addVertex(fusionUnmapped2RegexCatJob);
            graph.addEdge(mapspliceMultiThreadFusionJob, fusionUnmapped2RegexCatJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, ReadsToUnmappedSAMCLI.class, attempt.getId()).siteName(
                    siteName);
            File fusionUnmapped1ReadsToUnmappedSAMOutput = new File(fusionDir, "fusion_unmapped.1.sam");
            builder.addArgument(ReadsToUnmappedSAMCLI.INPUT, fusionUnmapped1RegexCatOutput.getAbsolutePath())
                    .addArgument(ReadsToUnmappedSAMCLI.OUTPUT,
                            fusionUnmapped1ReadsToUnmappedSAMOutput.getAbsolutePath());
            CondorJob fusionUnmapped1ReadsToUnmappedSAMJob = builder.build();
            logger.info(fusionUnmapped1ReadsToUnmappedSAMJob.toString());
            graph.addVertex(fusionUnmapped1ReadsToUnmappedSAMJob);
            graph.addEdge(fusionUnmapped1RegexCatJob, fusionUnmapped1ReadsToUnmappedSAMJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, ReadsToUnmappedSAMCLI.class, attempt.getId()).siteName(
                    siteName);
            File fusionUnmapped2ReadsToUnmappedSAMOutput = new File(fusionDir, "fusion_unmapped.2.sam");
            builder.addArgument(ReadsToUnmappedSAMCLI.INPUT, fusionUnmapped2RegexCatOutput.getAbsolutePath())
                    .addArgument(ReadsToUnmappedSAMCLI.OUTPUT,
                            fusionUnmapped2ReadsToUnmappedSAMOutput.getAbsolutePath());
            CondorJob fusionUnmapped2ReadsToUnmappedSAMJob = builder.build();
            logger.info(fusionUnmapped2ReadsToUnmappedSAMJob.toString());
            graph.addVertex(fusionUnmapped2ReadsToUnmappedSAMJob);
            graph.addEdge(fusionUnmapped2RegexCatJob, fusionUnmapped2ReadsToUnmappedSAMJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, AlignmentHandlerMultiCLI.class, attempt.getId())
                    .siteName(siteName).numberOfProcessors(8);
            filteredAlignmentBase = new File(fusionDir, "_filtered_fusion_alignments");
            builder.addArgument(AlignmentHandlerMultiCLI.ADDSOFTCLIP, 1)
                    .addArgument(AlignmentHandlerMultiCLI.AVERAGEFRAGMENTLENGTH, 225)
                    .addArgument(AlignmentHandlerMultiCLI.BOUNDARY, 36)
                    .addArgument(AlignmentHandlerMultiCLI.CHROMOSOMESIZEFILE, chromosomeSizesOutput.getAbsolutePath())
                    .addArgument(AlignmentHandlerMultiCLI.ENCOMPASSINGFUSIONREGIONEXTENSION, 50000)
                    .addArgument(AlignmentHandlerMultiCLI.FILTEREDALIGNMENTAPPEND, "_filtered_fusion_alignments")
                    .addArgument(AlignmentHandlerMultiCLI.FILTEREDALIGNMENTBASE,
                            filteredAlignmentBase.getAbsolutePath())
                    .addArgument(AlignmentHandlerMultiCLI.FILTERFLAG, 12 + 32 + 128 + 1024)
                    .addArgument(AlignmentHandlerMultiCLI.FRAGMENTLENGTH, 400)
                    .addArgument(AlignmentHandlerMultiCLI.FRAGMENTLENGTHSD, 100)
                    .addArgument(AlignmentHandlerMultiCLI.INPUT, sortOutput.getAbsolutePath())
                    .addArgument(AlignmentHandlerMultiCLI.INTRONDISTANCESD, 500)
                    .addArgument(AlignmentHandlerMultiCLI.MATEDISTANCESD, 100)
                    .addArgument(AlignmentHandlerMultiCLI.MAXIMUMANCHORDIFF, 50)
                    .addArgument(AlignmentHandlerMultiCLI.MAXIMUMDELETION, 6)
                    .addArgument(AlignmentHandlerMultiCLI.MAXIMUMHITS, 40)
                    .addArgument(AlignmentHandlerMultiCLI.MAXIMUMMATEDISTANCE, 50000)
                    .addArgument(AlignmentHandlerMultiCLI.READLENGTHPROPERTIES,
                            determineReadLengthOutput.getAbsolutePath())
                    .addArgument(AlignmentHandlerMultiCLI.MINIMUMANCHOR, 0)
                    .addArgument(AlignmentHandlerMultiCLI.MINIMUMCOVERAGE, 0)
                    .addArgument(AlignmentHandlerMultiCLI.MINIMUMENCOMPASSCOUNT, 1)
                    .addArgument(AlignmentHandlerMultiCLI.MINIMUMENTROPY, -0.0001)
                    .addArgument(AlignmentHandlerMultiCLI.MINIMUMINSERTION, 6)
                    .addArgument(AlignmentHandlerMultiCLI.MINIMUMJUNCTIONANCHOR, 10)
                    .addArgument(AlignmentHandlerMultiCLI.MINIMUMMISMATCH, 5)
                    .addArgument(AlignmentHandlerMultiCLI.THREADS, 8);
            alignmentHandlerMultiJob = builder.build();
            logger.info(alignmentHandlerMultiJob.toString());
            graph.addVertex(alignmentHandlerMultiJob);
            graph.addEdge(sortJob, alignmentHandlerMultiJob);
            graph.addEdge(determineReadLengthJob, alignmentHandlerMultiJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, SortCLI.class, attempt.getId()).siteName(siteName);
            sortOutput = new File(tmpDir, "combined_unmapped.sam");
            builder.addArgument(SortCLI.OUTPUT, sortOutput.getAbsolutePath())
                    .addArgument(SortCLI.BUFFERSIZE, 3500000)
                    .addArgument(SortCLI.KEY, "1,1")
                    .addArgument(SortCLI.TMPDIRECTORY, tmpDir.getAbsolutePath())
                    .addArgument(SortCLI.INPUT,
                            new File(remapDir, "_filtered_normal_alignments.unmapped").getAbsolutePath())
                    .addArgument(SortCLI.INPUT, fusionUnmapped1ReadsToUnmappedSAMOutput.getAbsolutePath())
                    .addArgument(SortCLI.INPUT, fusionUnmapped2ReadsToUnmappedSAMOutput.getAbsolutePath())
                    .addArgument(SortCLI.INPUT,
                            new File(fusionDir, "_filtered_fusion_alignments.unmapped").getAbsolutePath());
            sortJob = builder.build();
            logger.info(sortJob.toString());
            graph.addVertex(sortJob);
            graph.addEdge(alignmentHandlerMultiJob, sortJob);
            graph.addEdge(fusionUnmapped1ReadsToUnmappedSAMJob, sortJob);
            graph.addEdge(fusionUnmapped2ReadsToUnmappedSAMJob, sortJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, SetUnmappedBitFlagCLI.class, attempt.getId()).siteName(
                    siteName);
            File unmappedBitFlag = new File(tmpDir, "combined_unmapped_setbitflag.sam");
            builder.addArgument(SetUnmappedBitFlagCLI.UNMAPPEDSAM, sortOutput.getAbsolutePath()).addArgument(
                    SetUnmappedBitFlagCLI.UNMAPPEDSETBIT, unmappedBitFlag.getAbsolutePath());
            CondorJob setUnmappedBitFlagJob = builder.build();
            logger.info(setUnmappedBitFlagJob.toString());
            graph.addVertex(setUnmappedBitFlagJob);
            graph.addEdge(sortJob, setUnmappedBitFlagJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, CatCLI.class, attempt.getId()).siteName(siteName);
            File catOutput = new File(tmpDir, "final_alignments_headed.sam");
            builder.addArgument(CatCLI.OUTPUT, catOutput.getAbsolutePath())
                    .addArgument(CatCLI.FILES, chromosomeHeadOutput.getAbsolutePath())
                    .addArgument(CatCLI.FILES, pairedAlignments.getAbsolutePath())
                    .addArgument(CatCLI.FILES, singleAlignments.getAbsolutePath())
                    .addArgument(CatCLI.FILES, fusionPairedAlignments.getAbsolutePath())
                    .addArgument(CatCLI.FILES,
                            new File(remapDir, "_filtered_normal_alignments.unmapped").getAbsolutePath())
                    .addArgument(CatCLI.FILES, fusionUnmapped1ReadsToUnmappedSAMOutput.getAbsolutePath())
                    .addArgument(CatCLI.FILES, fusionUnmapped2ReadsToUnmappedSAMOutput.getAbsolutePath())
                    .addArgument(CatCLI.FILES,
                            new File(fusionDir, "_filtered_fusion_alignments.unmapped").getAbsolutePath())
                    .addArgument(CatCLI.FILES, unmappedBitFlag.getAbsolutePath());
            CondorJob finalAlignmentsHeadedCatJob = builder.build();
            logger.info(finalAlignmentsHeadedCatJob.toString());
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
            builder = WorkflowJobFactory.createJob(++count, SAMToolsViewCLI.class, attempt.getId()).siteName(siteName);
            File samtoolsViewOutput = new File(outputDirectory, fastqLaneRootName + ".bam");
            builder.addArgument(SAMToolsViewCLI.OUTPUT, samtoolsViewOutput.getAbsolutePath())
                    .addArgument(SAMToolsViewCLI.SAMINPUTFORMAT).addArgument(SAMToolsViewCLI.BAMFORMAT)
                    .addArgument(SAMToolsViewCLI.INPUT, catOutput.getAbsolutePath());
            CondorJob samtoolsViewJob = builder.build();
            logger.info(samtoolsViewJob.toString());
            graph.addVertex(samtoolsViewJob);
            graph.addEdge(finalAlignmentsHeadedCatJob, samtoolsViewJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, PicardAddOrReplaceReadGroupsCLI.class, attempt.getId(),
                    sample.getId()).siteName(siteName);
            File fixRGOutput = new File(outputDirectory, samtoolsViewOutput.getName().replace(".bam", ".fixed-rg.bam"));
            builder.addArgument(PicardAddOrReplaceReadGroupsCLI.INPUT, samtoolsViewOutput.getAbsolutePath())
                    .addArgument(PicardAddOrReplaceReadGroupsCLI.OUTPUT, fixRGOutput.getAbsolutePath())
                    .addArgument(PicardAddOrReplaceReadGroupsCLI.SORTORDER,
                            PicardSortOrderType.COORDINATE.toString().toLowerCase())
                    .addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPID, sample.getId().toString())
                    .addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPLIBRARY, sampleName)
                    .addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPPLATFORM, readGroupPlatform)
                    .addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPPLATFORMUNIT, readGroupPlatformUnit)
                    .addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPSAMPLENAME, sampleName);
            CondorJob picardAddOrReplaceReadGroupsJob = builder.build();
            logger.info(picardAddOrReplaceReadGroupsJob.toString());
            graph.addVertex(picardAddOrReplaceReadGroupsJob);
            graph.addEdge(samtoolsViewJob, picardAddOrReplaceReadGroupsJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, SAMToolsSortCLI.class, attempt.getId(), sample.getId())
                    .siteName(siteName);
            File samtoolsSortOut = new File(outputDirectory, fixRGOutput.getName().replace(".bam", ".sorted.bam"));
            builder.addArgument(SAMToolsSortCLI.INPUT, fixRGOutput.getAbsolutePath()).addArgument(
                    SAMToolsSortCLI.OUTPUT, samtoolsSortOut.getAbsolutePath());
            CondorJob samtoolsSortJob = builder.build();
            logger.info(samtoolsSortJob.toString());
            graph.addVertex(samtoolsSortJob);
            graph.addEdge(picardAddOrReplaceReadGroupsJob, samtoolsSortJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, SAMToolsIndexCLI.class, attempt.getId(), sample.getId())
                    .siteName(siteName);
            File samtoolsIndexOutput = new File(outputDirectory, samtoolsSortOut.getName().replace(".bam", ".bai"));
            builder.addArgument(SAMToolsIndexCLI.INPUT, samtoolsSortOut.getAbsolutePath()).addArgument(
                    SAMToolsIndexCLI.OUTPUT, samtoolsIndexOutput.getAbsolutePath());
            CondorJob samtoolsIndexJob = builder.build();
            logger.info(samtoolsIndexJob.toString());
            graph.addVertex(samtoolsIndexJob);
            graph.addEdge(samtoolsSortJob, samtoolsIndexJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, SAMToolsFlagstatCLI.class, attempt.getId(), sample.getId())
                    .siteName(siteName);
            File samtoolsFlagstatOut = new File(outputDirectory, samtoolsSortOut.getName().replace(".bam", ".flagstat"));
            builder.addArgument(SAMToolsFlagstatCLI.INPUT, samtoolsSortOut.getAbsolutePath()).addArgument(
                    SAMToolsFlagstatCLI.OUTPUT, samtoolsFlagstatOut.getAbsolutePath());
            CondorJob samtoolsFlagstatJob = builder.build();
            logger.info(samtoolsFlagstatJob.toString());
            graph.addVertex(samtoolsFlagstatJob);
            graph.addEdge(samtoolsIndexJob, samtoolsFlagstatJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, UBUSamJunctionCLI.class, attempt.getId(), sample.getId())
                    .siteName(siteName);
            File ubuSamJunctionOut = new File(outputDirectory, samtoolsSortOut.getName().replace(".bam",
                    ".junction_quantification.txt"));
            builder.addArgument(UBUSamJunctionCLI.JUNCTIONS, junctions)
                    .addArgument(UBUSamJunctionCLI.INPUT, samtoolsSortOut.getAbsolutePath())
                    .addArgument(UBUSamJunctionCLI.OUTPUT, ubuSamJunctionOut.getAbsolutePath());
            CondorJob ubuSamJunctionJob = builder.build();
            logger.info(ubuSamJunctionJob.toString());
            graph.addVertex(ubuSamJunctionJob);
            graph.addEdge(samtoolsIndexJob, ubuSamJunctionJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, CoverageBedCLI.class, attempt.getId(), sample.getId())
                    .siteName(siteName);
            File coverageBedOut = new File(outputDirectory, samtoolsSortOut.getName().replace(".bam",
                    ".coverageBedOut.txt"));
            builder.addArgument(CoverageBedCLI.INPUT, samtoolsSortOut.getAbsolutePath())
                    .addArgument(CoverageBedCLI.BED, compositeExons).addArgument(CoverageBedCLI.SPLITBED)
                    .addArgument(CoverageBedCLI.OUTPUT, coverageBedOut.getAbsolutePath());
            CondorJob coverageBedJob = builder.build();
            logger.info(coverageBedJob.toString());
            graph.addVertex(coverageBedJob);
            graph.addEdge(samtoolsSortJob, coverageBedJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, NormBedExonQuantCLI.class, attempt.getId(), sample.getId())
                    .siteName(siteName);
            File normBedExonQuantOut = new File(outputDirectory, coverageBedOut.getName().replace(
                    ".coverageBedOut.txt", ".normBedExonQuantOut.txt"));
            builder.addArgument(NormBedExonQuantCLI.INFILE, coverageBedOut.getAbsolutePath())
                    .addArgument(NormBedExonQuantCLI.COMPOSITEBED, compositeExons)
                    .addArgument(NormBedExonQuantCLI.OUTFILE, normBedExonQuantOut.getAbsolutePath());
            CondorJob normBedExonQuantJob = builder.build();
            logger.info(normBedExonQuantJob.toString());
            graph.addVertex(normBedExonQuantJob);
            graph.addEdge(coverageBedJob, normBedExonQuantJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, SortByReferenceAndNameCLI.class, attempt.getId(),
                    sample.getId()).siteName(siteName);
            File sortBAMByReferenceAndNameOut = new File(outputDirectory, samtoolsSortOut.getName().replace(".bam",
                    ".refAndNameSort.bam"));
            builder.addArgument(SortByReferenceAndNameCLI.INPUT, samtoolsSortOut.getAbsolutePath()).addArgument(
                    SortByReferenceAndNameCLI.OUTPUT, sortBAMByReferenceAndNameOut.getAbsolutePath());
            CondorJob sortBAMByReferenceAndNameJob = builder.build();
            logger.info(sortBAMByReferenceAndNameJob.toString());
            graph.addVertex(sortBAMByReferenceAndNameJob);
            graph.addEdge(samtoolsIndexJob, sortBAMByReferenceAndNameJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, UBUSamTranslateCLI.class, attempt.getId(), sample.getId())
                    .siteName(siteName);
            File ubuSamTranslateOut = new File(outputDirectory, sortBAMByReferenceAndNameOut.getName().replace(".bam",
                    ".transcriptomeAlignments.bam"));
            builder.addArgument(UBUSamTranslateCLI.BED, bed)
                    .addArgument(UBUSamTranslateCLI.INPUT, sortBAMByReferenceAndNameOut.getAbsolutePath())
                    .addArgument(UBUSamTranslateCLI.OUTPUT, ubuSamTranslateOut.getAbsolutePath())
                    .addArgument(UBUSamTranslateCLI.ORDER, getWorkflowBeanService().getAttributes().get("order"))
                    .addArgument(UBUSamTranslateCLI.XGTAGS).addArgument(UBUSamTranslateCLI.REVERSE);
            CondorJob ubuSamTranslateJob = builder.build();
            logger.info(ubuSamTranslateJob.toString());
            graph.addVertex(ubuSamTranslateJob);
            graph.addEdge(sortBAMByReferenceAndNameJob, ubuSamTranslateJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, UBUSamFilterCLI.class, attempt.getId(), sample.getId())
                    .siteName(siteName);
            File ubuSamFilterOut = new File(outputDirectory, ubuSamTranslateOut.getName().replace(".bam",
                    ".filtered.bam"));
            builder.addArgument(UBUSamFilterCLI.INPUT, ubuSamTranslateOut.getAbsolutePath())
                    .addArgument(UBUSamFilterCLI.OUTPUT, ubuSamFilterOut.getAbsolutePath())
                    .addArgument(UBUSamFilterCLI.STRIPINDELS).addArgument(UBUSamFilterCLI.MAXINSERT, "10000")
                    .addArgument(UBUSamFilterCLI.MAPQ, "1");
            CondorJob ubuSamFilterJob = builder.build();
            logger.info(ubuSamFilterJob.toString());
            graph.addVertex(ubuSamFilterJob);
            graph.addEdge(ubuSamTranslateJob, ubuSamFilterJob);

            // new job
            builder = WorkflowJobFactory
                    .createJob(++count, RSEMCalculateExpressionCLI.class, attempt.getId(), sample.getId())
                    .siteName(siteName).numberOfProcessors(4);
            File outputPrefix = new File(outputDirectory, "rsem");
            builder.addArgument(RSEMCalculateExpressionCLI.PAIREDEND).addArgument(RSEMCalculateExpressionCLI.BAM)
                    .addArgument(RSEMCalculateExpressionCLI.ESTIMATERSPD)
                    .addArgument(RSEMCalculateExpressionCLI.THREADS, "4")
                    .addArgument(RSEMCalculateExpressionCLI.BAMFILE, ubuSamFilterOut.getAbsolutePath())
                    .addArgument(RSEMCalculateExpressionCLI.REFERENCESEQUENCE, referenceSequencePrefix)
                    .addArgument(RSEMCalculateExpressionCLI.OUTPUT, outputPrefix.getAbsolutePath());
            CondorJob rsemJob = builder.build();
            logger.info(rsemJob.toString());
            graph.addVertex(rsemJob);
            graph.addEdge(ubuSamFilterJob, rsemJob);

            // new job
            builder = WorkflowJobFactory
                    .createJob(++count, StripTrailingTabsCLI.class, attempt.getId(), sample.getId()).siteName(siteName);
            File isoformResults = new File(outputDirectory, outputPrefix.getName() + ".isoforms.results");
            File isoformResultsStripped = new File(outputDirectory, outputPrefix.getName()
                    + ".isoforms.stripped.results");
            builder.addArgument(StripTrailingTabsCLI.INPUT, isoformResults.getAbsolutePath()).addArgument(
                    StripTrailingTabsCLI.OUTPUT, isoformResultsStripped.getAbsolutePath());
            builder.postScript(String.format("/bin/mv %s %s", isoformResultsStripped.getAbsolutePath(),
                    isoformResults.getAbsolutePath()));
            CondorJob rsemISOFormsResultStripTabJob = builder.build();
            logger.info(rsemISOFormsResultStripTabJob.toString());
            graph.addVertex(rsemISOFormsResultStripTabJob);
            graph.addEdge(rsemJob, rsemISOFormsResultStripTabJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, PruneISOFormsFromGeneQuantFileCLI.class, attempt.getId(),
                    sample.getId()).siteName(siteName);
            File geneResultsFile = new File(outputDirectory, outputPrefix.getName() + ".genes.results");
            File origGeneResultsFile = new File(outputDirectory, "orig." + outputPrefix.getName() + ".genes.results");
            builder.addArgument(PruneISOFormsFromGeneQuantFileCLI.GENERESULTS, geneResultsFile.getAbsolutePath())
                    .addArgument(PruneISOFormsFromGeneQuantFileCLI.ORIGGENERESULTS,
                            origGeneResultsFile.getAbsolutePath());
            CondorJob pruneISOFormsFromGeneQuantFileJob = builder.build();
            logger.info(pruneISOFormsFromGeneQuantFileJob.toString());
            graph.addVertex(pruneISOFormsFromGeneQuantFileJob);
            graph.addEdge(rsemISOFormsResultStripTabJob, pruneISOFormsFromGeneQuantFileJob);

            // new job
            builder = WorkflowJobFactory
                    .createJob(++count, NormalizeQuartileCLI.class, attempt.getId(), sample.getId()).siteName(siteName);
            File normalizeGeneQuantOut = new File(outputDirectory, outputPrefix.getName() + ".genes.normalized_results");
            builder.addArgument(NormalizeQuartileCLI.COLUMN, "2").addArgument(NormalizeQuartileCLI.QUANTILE, "75")
                    .addArgument(NormalizeQuartileCLI.TARGET, "1000")
                    .addArgument(NormalizeQuartileCLI.INPUT, geneResultsFile.getAbsolutePath())
                    .addArgument(NormalizeQuartileCLI.OUTPUT, normalizeGeneQuantOut.getAbsolutePath());
            CondorJob normalizeGeneQuantJob = builder.build();
            logger.info(normalizeGeneQuantJob.toString());
            graph.addVertex(normalizeGeneQuantJob);
            graph.addEdge(pruneISOFormsFromGeneQuantFileJob, normalizeGeneQuantJob);

            // new job
            builder = WorkflowJobFactory
                    .createJob(++count, NormalizeQuartileCLI.class, attempt.getId(), sample.getId()).siteName(siteName);
            File normalizeISOFormQuantOut = new File(outputDirectory, outputPrefix.getName()
                    + ".isoforms.normalized_results");
            builder.addArgument(NormalizeQuartileCLI.COLUMN, "2").addArgument(NormalizeQuartileCLI.QUANTILE, "75")
                    .addArgument(NormalizeQuartileCLI.TARGET, "300")
                    .addArgument(NormalizeQuartileCLI.INPUT, isoformResults.getAbsolutePath())
                    .addArgument(NormalizeQuartileCLI.OUTPUT, normalizeISOFormQuantOut.getAbsolutePath());
            CondorJob normalizeISOFormQuantJob = builder.build();
            logger.info(normalizeISOFormQuantJob.toString());
            graph.addVertex(normalizeISOFormQuantJob);
            graph.addEdge(pruneISOFormsFromGeneQuantFileJob, normalizeISOFormQuantJob);

            // new job
            builder = WorkflowJobFactory.createJob(++count, RemoveCLI.class, attempt.getId(), sample.getId()).siteName(
                    siteName);
            builder.addArgument(RemoveCLI.FILE, gunzippedFastqR1.getAbsolutePath())
                    .addArgument(RemoveCLI.FILE, gunzippedFastqR2.getAbsolutePath())
                    .addArgument(RemoveCLI.FILE, tmpDir.getAbsolutePath())
                    .addArgument(RemoveCLI.FILE, fastqFormatterR1Out.getAbsolutePath())
                    .addArgument(RemoveCLI.FILE, fastqFormatterR2Out.getAbsolutePath())
                    .addArgument(RemoveCLI.FILE, samtoolsViewOutput.getAbsolutePath())
                    .addArgument(RemoveCLI.FILE, fixRGOutput.getAbsolutePath());
            CondorJob removeJob = builder.build();
            logger.info(removeJob.toString());
            graph.addVertex(removeJob);
            graph.addEdge(normalizeISOFormQuantJob, removeJob);
            graph.addEdge(normalizeGeneQuantJob, removeJob);

        }

        return graph;
    }
}
