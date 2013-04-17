package edu.unc.mapseq.messaging.rnaseq;

import java.util.Date;

import javax.jms.Connection;
import javax.jms.DeliveryMode;
import javax.jms.Destination;
import javax.jms.JMSException;
import javax.jms.MessageProducer;
import javax.jms.Session;

import org.apache.activemq.ActiveMQConnectionFactory;
import org.apache.commons.lang.time.DateFormatUtils;
import org.junit.Test;

public class RNASeqMessageTest {

    @Test
    public void testQueue() {

        ActiveMQConnectionFactory connectionFactory = new ActiveMQConnectionFactory(String.format("nio://%s:61616",
                "gnet641.its.unc.edu"));

        Connection connection = null;
        Session session = null;
        try {
            connection = connectionFactory.createConnection();
            session = connection.createSession(false, Session.AUTO_ACKNOWLEDGE);
            Destination destination = session.createQueue("queue/rnaseq");
            MessageProducer producer = session.createProducer(destination);
            producer.setDeliveryMode(DeliveryMode.NON_PERSISTENT);
            String format = "{\"account_name\":\"rc_lbg.svc\",\"entities\":[{\"entity_type\":\"HTSFSample\",\"guid\":\"%1$d\"},{\"entity_type\":\"WorkflowRun\",\"name\":\"test-%2$s-%1$d\"}]}";
            producer.send(session.createTextMessage(String.format(format, 314616,
                    DateFormatUtils.ISO_DATE_FORMAT.format(new Date()))));
        } catch (JMSException e) {
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

}
