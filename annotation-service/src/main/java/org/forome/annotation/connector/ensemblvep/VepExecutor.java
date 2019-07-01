package org.forome.annotation.connector.ensemblvep;

import com.jcraft.jsch.*;
import org.forome.annotation.config.connector.EnsemblVepConfigConnector;
import org.forome.annotation.config.sshtunnel.SshTunnelConfig;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;

public class VepExecutor implements Closeable {

    private final static Logger log = LoggerFactory.getLogger(VepExecutor.class);

    private final EnsemblVepConfigConnector ensemblVepConfigConnector;

    private JSch jsch;
    private Session sshSession;

    public VepExecutor(EnsemblVepConfigConnector ensemblVepConfigConnector) throws JSchException {
        this.ensemblVepConfigConnector = ensemblVepConfigConnector;
        connect();
    }

    public synchronized void connect() throws JSchException {
        SshTunnelConfig sshTunnelConfig = ensemblVepConfigConnector.sshTunnelConfig;
        if (sshTunnelConfig != null) {

            jsch = new JSch();
            jsch.addIdentity(sshTunnelConfig.key);

            sshSession = jsch.getSession(sshTunnelConfig.user, sshTunnelConfig.host, sshTunnelConfig.port);

            java.util.Properties config = new java.util.Properties();
            config.put("StrictHostKeyChecking", "no");
            config.put("Compression", "yes");
            config.put("ConnectionAttempts", "2");
            sshSession.setConfig(config);

            sshSession.connect();
        }
    }

    private synchronized void disconnect() {
        if (sshSession != null) {
            sshSession.disconnect();
            sshSession = null;
        }
        jsch = null;
    }

    public String exec(String command, InputStream stdin) throws Exception {

        Channel channel = sshSession.openChannel("exec");
        ((ChannelExec)channel).setCommand(command);
        channel.setInputStream(stdin);
        ((ChannelExec)channel).setErrStream(System.err);

        InputStream input = channel.getInputStream();
        channel.connect();

        try{
            InputStreamReader inputReader = new InputStreamReader(input);
            BufferedReader bufferedReader = new BufferedReader(inputReader);
            String line = null;

            while((line = bufferedReader.readLine()) != null){
                System.out.println(line);
            }
            bufferedReader.close();
            inputReader.close();
        }catch(IOException ex){
            ex.printStackTrace();
        }

        channel.disconnect();

        return null;
    }

    @Override
    public void close() {
        disconnect();
    }
}
