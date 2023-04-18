# Functions for logging the solver output

# Function for printing a message to stdout and a log file
function report(message, returnTime = 0)
    # Get the current time
    now = Dates.now();
    # Prepend the current time to the message
    message = "$(Dates.format(now, "yyyy-mm-dd HH:MM:SS.sss"))   $(message)";
    # Print the message to stdout
    println(message)
    # Write the message to the log file
    file = open("simulation.log", "a");
    println(file, message);
    close(file)
    if (returnTime == 1)
        return now
    end
end