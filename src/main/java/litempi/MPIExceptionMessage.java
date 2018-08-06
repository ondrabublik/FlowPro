package litempi;

/**
 *
 * @author ales
 */
public class MPIExceptionMessage extends MPIMessage {
    private final Throwable ex;
    
    public MPIExceptionMessage(Throwable ex) {
        super(Integer.MIN_VALUE);
        this.ex = ex;
    }
    
    public Throwable getException() {
        return ex;
    }
}
