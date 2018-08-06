package litempi;

import java.io.Serializable;

public class MPIMessage implements Serializable {
    public final int tag;
    private final Object obj;
    
    public MPIMessage(int tag, Object obj) {
        this.tag = tag;
        this.obj = obj;
    }
    
    public MPIMessage(int tag) {
        this.tag = tag;
        obj = null;
    }
    
    public Object getData() {
        return obj;
    }
}
