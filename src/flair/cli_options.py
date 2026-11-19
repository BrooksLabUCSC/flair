"""Options whose values may be comma separated.  Such an option takes its values
comma separated, by repeating the option, or both, so a long list needn't be one
argument of the command line.  An option whose values may contain a comma, such as
--region, where a position may be written 1,000, takes them only by repetition."""

def option_values(values):
    "the values of a list option, however they were given; empty when it was not"
    return [value.strip() for group in (values or ()) for value in group.split(',') if value.strip() != '']
