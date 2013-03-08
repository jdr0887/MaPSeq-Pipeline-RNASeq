package edu.unc.mapseq.messaging.rnaseq;

import javax.jms.Connection;
import javax.jms.ConnectionFactory;
import javax.jms.Destination;
import javax.jms.MessageConsumer;
import javax.jms.Session;

import org.osgi.framework.BundleActivator;
import org.osgi.framework.BundleContext;
import org.osgi.framework.ServiceReference;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.pipeline.PipelineBeanService;
import edu.unc.mapseq.pipeline.rnaseq.RNASeqPipelineBeanService;

public class Activator implements BundleActivator {

    private final Logger logger = LoggerFactory.getLogger(Activator.class);

    private Connection connection;

    private Session session;

    private Destination destination;

    @Override
    public void start(BundleContext context) throws Exception {
        logger.debug("ENTERING start(BundleContext)");

        ServiceReference reference = context.getServiceReference(ConnectionFactory.class.getName());
        ConnectionFactory connectionFactory = (ConnectionFactory) context.getService(reference);

        RNASeqPipelineBeanService pipelineBeanService = null;
        ServiceReference[] references = context.getServiceReferences(PipelineBeanService.class.getName(), null);
        if (references != null) {
            for (ServiceReference ref : references) {
                Object o = context.getService(ref);
                if (o instanceof RNASeqPipelineBeanService) {
                    pipelineBeanService = (RNASeqPipelineBeanService) o;
                    break;
                }
            }
        }

        if (pipelineBeanService != null) {
            this.connection = connectionFactory.createConnection();
            this.session = connection.createSession(false, Session.AUTO_ACKNOWLEDGE);
            this.destination = this.session.createQueue("queue/rnaseq");
            MessageConsumer consumer = session.createConsumer(destination);
            RNASeqMessageListener messageListener = new RNASeqMessageListener();
            messageListener.setPipelineBeanService(pipelineBeanService);
            consumer.setMessageListener(messageListener);
            this.connection.start();
        }

    }

    @Override
    public void stop(BundleContext context) throws Exception {
        logger.debug("ENTERING stop(BundleContext)");
        if (this.session != null) {
            this.session.close();
        }
        if (this.connection != null) {
            this.connection.stop();
            this.connection.close();
        }
    }

}
