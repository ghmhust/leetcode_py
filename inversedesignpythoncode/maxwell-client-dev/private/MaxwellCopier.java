import java.io.*;
import java.lang.Math;

public class MaxwellCopier{
    // This java class allows for simultaneous communication across HTTP 
    // for Maxwell.
    //
    // TODO: Make this permanent or something
    // private static byte[] buf = new byte[4 * 1024];
    public int total_bytes_transferred;
    private byte[] buf;

    public MaxwellCopier() { 
        total_bytes_transferred = 0; 
        buf = new byte[4 * 1024];
        }

    public boolean copy(InputStream input, OutputStream output) throws Exception {
        // int len = input.read(buf, 0, Math.max(1, Math.min(input.available(), 1 * 1024)));
        int len = input.read(buf);
        if (len == -1) {
            return false; // No longer running.
        } else {
            total_bytes_transferred += len;
            output.write(buf, 0, len);
            return true; // Still running.
        }
    }
}
