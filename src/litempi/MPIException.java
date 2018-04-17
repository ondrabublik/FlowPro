package litempi;

/**
 *
 * @author ales
 */
public class MPIException extends Exception {

    public MPIException(String message) {
        super(message);
    }

    public MPIException(String message, Throwable throwable) {
        super(message, throwable);
    }
}
