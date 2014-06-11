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

import edu.unc.mapseq.dao.MaPSeqDAOException;
import edu.unc.mapseq.dao.model.Account;
import edu.unc.mapseq.dao.model.HTSFSample;
import edu.unc.mapseq.dao.model.SequencerRun;
import edu.unc.mapseq.dao.ws.WSDAOManager;

public class RunWorkflow implements Runnable {

    private final Logger logger = LoggerFactory.getLogger(RunWorkflow.class);

    private final static HelpFormatter helpFormatter = new HelpFormatter();

    private final static Options cliOptions = new Options();

    private final WSDAOManager daoMgr = WSDAOManager.getInstance();

    private String workflowRunName;

    private Long sequencerRunId;

    private Long htsfSampleId;

    private String host = "gnet641.its.unc.edu";

    public RunWorkflow() {
        super();
    }

    @Override
    public void run() {
        logger.info("ENTERING run()");
        Account account = null;
        try {
            account = daoMgr.getMaPSeqDAOBean().getAccountDAO().findByName(System.getProperty("user.name"));
        } catch (MaPSeqDAOException e) {
            e.printStackTrace();
        }

        if (account == null) {
            System.out.println("Must register account first");
            return;
        }

        if (this.sequencerRunId == null && this.htsfSampleId == null) {
            System.err.println("Both sequencerRunId and htsfSampeId can't be null");
            return;
        }

        SequencerRun sequencerRun = null;
        Set<HTSFSample> htsfSampleSet = new HashSet<HTSFSample>();

        if (sequencerRunId != null && htsfSampleId == null) {
            try {
                sequencerRun = daoMgr.getMaPSeqDAOBean().getSequencerRunDAO().findById(this.sequencerRunId);
            } catch (MaPSeqDAOException e) {
                logger.error("Problem finding SequencerRun", e);
            }
        } else if (sequencerRunId == null && htsfSampleId != null) {
            try {
                HTSFSample htsfSample = daoMgr.getMaPSeqDAOBean().getHTSFSampleDAO().findById(this.htsfSampleId);
                htsfSampleSet.add(htsfSample);
            } catch (MaPSeqDAOException e) {
                logger.error("Problem finding HTSFSample", e);
            }
        }

        if (sequencerRun == null && htsfSampleSet.size() == 0) {
            System.err.println("SequencerRun & Set<HTSFSample> are both null or empty");
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
            generator.writeStringField("accountName", System.getProperty("user.name"));
            generator.writeArrayFieldStart("entities");

            if (htsfSampleSet != null) {
                for (HTSFSample sample : htsfSampleSet) {
                    generator.writeStartObject();
                    generator.writeStringField("entityType", "HTSFSample");
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

    public Long getSequencerRunId() {
        return sequencerRunId;
    }

    public void setSequencerRunId(Long sequencerRunId) {
        this.sequencerRunId = sequencerRunId;
    }

    public Long getHtsfSampleId() {
        return htsfSampleId;
    }

    public void setHtsfSampleId(Long htsfSampleId) {
        this.htsfSampleId = htsfSampleId;
    }

    @SuppressWarnings("static-access")
    public static void main(String[] args) {
        cliOptions.addOption(OptionBuilder.withArgName("htsfSampleId").hasArg()
                .withDescription("HTSFSample identifier").withLongOpt("htsfSampleId").create());
        cliOptions.addOption(OptionBuilder.withArgName("sequencerRunId").hasArg()
                .withDescription("SequencerRun identifier").withLongOpt("sequencerRunId").create());
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

            if (commandLine.hasOption("sequencerRunId")) {
                Long sequencerRunId = Long.valueOf(commandLine.getOptionValue("sequencerRunId"));
                main.setSequencerRunId(sequencerRunId);
            }

            if (commandLine.hasOption("htsfSampleId")) {
                Long htsfSampleId = Long.valueOf(commandLine.getOptionValue("htsfSampleId"));
                main.setHtsfSampleId(htsfSampleId);
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
