package edu.unc.mapseq.rnaseq;

import java.io.IOException;
import java.io.StringWriter;
import java.util.HashSet;
import java.util.Set;

import javax.jms.Connection;
import javax.jms.DeliveryMode;
import javax.jms.Destination;
import javax.jms.JMSException;
import javax.jms.MessageProducer;
import javax.jms.Session;

import org.apache.activemq.ActiveMQConnectionFactory;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.fasterxml.jackson.core.JsonFactory;
import com.fasterxml.jackson.core.JsonGenerator;

import edu.unc.mapseq.dao.FlowcellDAO;
import edu.unc.mapseq.dao.MaPSeqDAOBean;
import edu.unc.mapseq.dao.MaPSeqDAOException;
import edu.unc.mapseq.dao.SampleDAO;
import edu.unc.mapseq.dao.model.Flowcell;
import edu.unc.mapseq.dao.model.Sample;
import edu.unc.mapseq.dao.ws.WSDAOManager;

public class RunWorkflow implements Runnable {

    private final Logger logger = LoggerFactory.getLogger(RunWorkflow.class);

    private final static HelpFormatter helpFormatter = new HelpFormatter();

    private final static Options cliOptions = new Options();

    private final WSDAOManager daoMgr = WSDAOManager.getInstance();

    private String workflowRunName;

    private Long flowcellId;

    private Long sampleId;

    private String host = "gnet641.its.unc.edu";

    public RunWorkflow() {
        super();
    }

    @Override
    public void run() {
        logger.info("ENTERING run()");

        if (this.flowcellId == null && this.sampleId == null) {
            System.err.println("Both flowcellId and sampleId can't be null");
            return;
        }

        Flowcell flowcell = null;
        Set<Sample> sampleSet = new HashSet<Sample>();

        MaPSeqDAOBean mapseqDAOBean = daoMgr.getMaPSeqDAOBean();
        FlowcellDAO flowcellDAO = mapseqDAOBean.getFlowcellDAO();
        SampleDAO sampleDAO = mapseqDAOBean.getSampleDAO();

        if (flowcellId != null && flowcellId == null) {
            try {
                flowcell = flowcellDAO.findById(this.flowcellId);
            } catch (MaPSeqDAOException e) {
                logger.error("Problem finding SequencerRun", e);
            }
        } else if (flowcellId == null && sampleId != null) {
            try {
                Sample sample = sampleDAO.findById(this.sampleId);
                sampleSet.add(sample);
            } catch (MaPSeqDAOException e) {
                logger.error("Problem finding HTSFSample", e);
            }
        }

        if (flowcell == null && sampleSet.size() == 0) {
            System.err.println("Flowcell & Set<Sample> are both null or empty");
            return;
        }

        ActiveMQConnectionFactory connectionFactory = new ActiveMQConnectionFactory(String.format("nio://%s:61616",
                this.host));

        Connection connection = null;
        Session session = null;
        try {
            connection = connectionFactory.createConnection();
            session = connection.createSession(false, Session.AUTO_ACKNOWLEDGE);
            Destination destination = session.createQueue("queue/rnaseq");
            MessageProducer producer = session.createProducer(destination);
            producer.setDeliveryMode(DeliveryMode.PERSISTENT);

            StringWriter sw = new StringWriter();

            JsonGenerator generator = new JsonFactory().createGenerator(sw);

            generator.writeStartObject();
            generator.writeArrayFieldStart("entities");

            if (sampleSet != null) {
                for (Sample sample : sampleSet) {
                    generator.writeStartObject();
                    generator.writeStringField("entityType", "Sample");
                    generator.writeStringField("guid", sample.getId().toString());
                    generator.writeEndObject();
                }
            }

            generator.writeStartObject();
            generator.writeStringField("entityType", "WorkflowRun");
            generator.writeStringField("name", workflowRunName);
            generator.writeEndObject();

            generator.writeEndArray();
            generator.writeEndObject();

            generator.flush();
            generator.close();

            sw.flush();
            sw.close();

            producer.send(session.createTextMessage(sw.toString()));

        } catch (JMSException | IOException e) {
            e.printStackTrace();
        } finally {
            try {
                session.close();
                connection.close();
            } catch (JMSException e) {
                e.printStackTrace();
            }
        }

    }

    public String getHost() {
        return host;
    }

    public void setHost(String host) {
        this.host = host;
    }

    public String getWorkflowRunName() {
        return workflowRunName;
    }

    public void setWorkflowRunName(String workflowRunName) {
        this.workflowRunName = workflowRunName;
    }

    public Long getFlowcellId() {
        return flowcellId;
    }

    public void setFlowcellId(Long flowcellId) {
        this.flowcellId = flowcellId;
    }

    public Long getSampleId() {
        return sampleId;
    }

    public void setSampleId(Long sampleId) {
        this.sampleId = sampleId;
    }

    @SuppressWarnings("static-access")
    public static void main(String[] args) {
        cliOptions.addOption(OptionBuilder.withArgName("sampleId").hasArg().withDescription("HTSFSample identifier")
                .withLongOpt("sampleId").create());
        cliOptions.addOption(OptionBuilder.withArgName("flowcellId").hasArg()
                .withDescription("SequencerRun identifier").withLongOpt("flowcellId").create());
        cliOptions.addOption(OptionBuilder.withArgName("workflowRunName").withLongOpt("workflowRunName").isRequired()
                .hasArg().create());
        cliOptions.addOption(OptionBuilder.withArgName("host").withLongOpt("host").hasArg().create());

        RunWorkflow main = new RunWorkflow();
        CommandLineParser commandLineParser = new GnuParser();
        try {
            CommandLine commandLine = commandLineParser.parse(cliOptions, args);
            if (commandLine.hasOption("?")) {
                helpFormatter.printHelp(main.getClass().getSimpleName(), cliOptions);
                return;
            }

            if (commandLine.hasOption("workflowRunName")) {
                String workflowRunName = commandLine.getOptionValue("workflowRunName");
                main.setWorkflowRunName(workflowRunName);
            }

            if (commandLine.hasOption("flowcellId")) {
                Long flowcellId = Long.valueOf(commandLine.getOptionValue("flowcellId"));
                main.setFlowcellId(flowcellId);
            }

            if (commandLine.hasOption("sampleId")) {
                Long sampleId = Long.valueOf(commandLine.getOptionValue("sampleId"));
                main.setSampleId(sampleId);
            }

            if (commandLine.hasOption("host")) {
                main.setHost(commandLine.getOptionValue("host"));
            }

            main.run();
        } catch (ParseException e) {
            System.err.println(("Parsing Failed: " + e.getMessage()));
            helpFormatter.printHelp(main.getClass().getSimpleName(), cliOptions);
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println(("Error: " + e.getMessage()));
            helpFormatter.printHelp(main.getClass().getSimpleName(), cliOptions);
        }

    }
}
