package litempi;

import java.io.Serializable;

public class MPIMessage implements Serializable {
    public final int tag, subTag;
    private final Object obj;
    
    public MPIMessage(int tag, Object obj) {
        this.tag = tag;
        this.subTag = 0;
        this.obj = obj;
    }
    
    public MPIMessage(int tag) {
        this.tag = tag;
        this.subTag = 0;
        obj = null;
    }
    
    public MPIMessage(int tag, int subTag, Object obj) {
        this.tag = tag;
        this.subTag = subTag;
        this.obj = obj;
    }
    
    public MPIMessage(int tag, int subTag) {
        this.tag = tag;
        this.subTag = subTag;
        obj = null;
    }
    
    public Object getData() {
        return obj;
    }
}
