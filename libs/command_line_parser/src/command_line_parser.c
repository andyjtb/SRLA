#include "command_line_parser.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>

/* Measure the spec list size */
static uint32_t CommandLineParser_GetNumSpecifications(
        const struct CommandLineParserSpecification* clps)
{
    uint32_t num_specs;

    assert(clps != NULL);

    /* Advance the pointer until it hits the end of the list, 0 */
    num_specs = 0;
    while ((clps->short_option != 0) || (clps->long_option != NULL)) {
        num_specs++;
        clps++;
    }

    return num_specs;
}

/* Check command line parser specs */
static CommandLineParserBool CommandLineParser_CheckSpecification(
        const struct CommandLineParserSpecification* clps)
{
    uint32_t spec_no;
    uint32_t num_specs;

    assert(clps != NULL);

    /* Get the number of specifications */
    num_specs = CommandLineParser_GetNumSpecifications(clps);

    for (spec_no = 0; spec_no < num_specs; spec_no++) {
        uint32_t j;
        for (j = 0; j < num_specs; j++) {
            if (j == spec_no) {
                continue;
            }
            /* It is invalid if there are any options with the same string */
            if (clps[j].short_option == clps[spec_no].short_option) {
                return COMMAND_LINE_PARSER_FALSE;
            } else if ((clps[j].long_option != NULL) && (clps[spec_no].long_option != NULL)) {
                if (strcmp(clps[j].long_option, clps[spec_no].long_option) == 0) {
                    return COMMAND_LINE_PARSER_FALSE;
                }
            }
        }
    }

    /* No problem */
    return COMMAND_LINE_PARSER_TRUE;
}

/* Print argument description */
void CommandLineParser_PrintDescription(const struct CommandLineParserSpecification* clps)
{
    uint32_t  spec_no;
    char      arg_option_attr[256];
    char      command_str[256];
    uint32_t  num_specs;

    /* Argument check */
    if (clps == NULL) {
        fprintf(stderr, "Pointer to command-line specification is NULL. \n");
        return;
    }

    /* Check the specs */
    if (CommandLineParser_CheckSpecification(clps) != COMMAND_LINE_PARSER_TRUE) {
        fprintf(stderr, "Warning: Command-line specification is invalid. (Unable to parse) \n");
    }

    /* Get the number of specifications */
    num_specs = CommandLineParser_GetNumSpecifications(clps);

    /* Display specifications in order */
    for (spec_no = 0; spec_no < num_specs; spec_no++) {
        const struct CommandLineParserSpecification* pspec = &clps[spec_no];
        /* Create an attribute string for the arguments */
        if (pspec->need_argument == COMMAND_LINE_PARSER_TRUE) {
            sprintf(arg_option_attr, "(needs argument)");
        } else {
            strcpy(arg_option_attr, "");
        }

        /* Create a command string */
        if (pspec->long_option != NULL) {
            sprintf(command_str, "  -%c, --%s", pspec->short_option, pspec->long_option);
        } else {
            sprintf(command_str, "  -%c", pspec->short_option);
        }

        /* Print everything with explanation */
        printf("%-20s %-18s  %s \n",
                command_str, arg_option_attr,
                (pspec->description != NULL) ? pspec->description : "");
    }
}

/* Get index from option name */
static CommandLineParserResult CommandLineParser_GetSpecificationIndex(
        const struct CommandLineParserSpecification* clps,
        const char* option_name, uint32_t* index)
{
    uint32_t spec_no;
    uint32_t num_specs;

    /* Argument check */
    if (clps == NULL || option_name == NULL || index == NULL) {
        return COMMAND_LINE_PARSER_RESULT_INVALID_ARGUMENT;
    }

    /* Get the number of specifications */
    num_specs = CommandLineParser_GetNumSpecifications(clps);

    /* Search from short options */
    if (strlen(option_name) == 1) {
        for (spec_no = 0; spec_no < num_specs; spec_no++) {
            if (option_name[0] == clps[spec_no].short_option) {
                *index = spec_no;
                return COMMAND_LINE_PARSER_RESULT_OK;
            }
        }
    }

    /* Search from long option */
    for (spec_no = 0; spec_no < num_specs; spec_no++) {
        if (strcmp(option_name, clps[spec_no].long_option) == 0) {
            *index = spec_no;
            return COMMAND_LINE_PARSER_RESULT_OK;
        }
    }

    /* Not Found */
    return COMMAND_LINE_PARSER_RESULT_UNKNOWN_OPTION;
}

/* Get whether the option was specified from the option name */
CommandLineParserBool CommandLineParser_GetOptionAcquired(
        const struct CommandLineParserSpecification* clps,
        const char* option_name)
{
    uint32_t spec_no;

    /* Get index */
    if (CommandLineParser_GetSpecificationIndex(clps, option_name, &spec_no) != COMMAND_LINE_PARSER_RESULT_OK) {
        return COMMAND_LINE_PARSER_FALSE;
    }

    return clps[spec_no].acquired;
}

/* Get the option arguments from the option name */
const char* CommandLineParser_GetArgumentString(
        const struct CommandLineParserSpecification* clps,
        const char* option_name)
{
    uint32_t spec_no;

    /* Get index */
    if (CommandLineParser_GetSpecificationIndex(clps, option_name, &spec_no) != COMMAND_LINE_PARSER_RESULT_OK) {
        return NULL;
    }

    return clps[spec_no].argument_string;
}

/* Parse arguments */
CommandLineParserResult CommandLineParser_ParseArguments(
        struct CommandLineParserSpecification* clps,
        int32_t argc, const char* const* argv,
        const char** other_string_array, uint32_t other_string_array_size)
{
    int32_t     count;
    uint32_t    spec_no;
    const char* arg_str;
    uint32_t    other_string_index;
    uint32_t    num_specs;

    /* Argument check */
    if (argv == NULL || clps == NULL) {
        return COMMAND_LINE_PARSER_RESULT_INVALID_ARGUMENT;
    }

    /* Get the number of specifications */
    num_specs = CommandLineParser_GetNumSpecifications(clps);

    /* Check command line specifications */
    if (CommandLineParser_CheckSpecification(clps) != COMMAND_LINE_PARSER_TRUE) {
        return COMMAND_LINE_PARSER_RESULT_INVALID_SPECIFICATION;
    }

    /* Set all options to unacquired state */
    for (spec_no = 0; spec_no < num_specs; spec_no++) {
        clps[spec_no].acquired = COMMAND_LINE_PARSER_FALSE;
    }

    /* argv[0] is the program name, so skip it */
    other_string_index = 0;
    for (count = 1; count < argc; count++) {
        /* Get the elements of a string array */
        arg_str = argv[count];
        /* Check the option string */
        if (strncmp(arg_str, "--", 2) == 0) {
            /* Long option */
            for (spec_no = 0; spec_no < num_specs; spec_no++) {
                uint32_t long_option_len;
                struct CommandLineParserSpecification* pspec = &clps[spec_no];
                /* Skip if long option string is NULL */
                if (pspec->long_option == NULL) {
                    continue;
                }
                long_option_len = (uint32_t)strlen(pspec->long_option);
                if (strncmp(&arg_str[2], pspec->long_option, long_option_len) == 0) {
                    /* Long options must be followed by a null terminator or an '=' to specify the option */
                    if (arg_str[2 + long_option_len] == '\0') {
                        /* An option that has already been acquired was specified */
                        if (pspec->acquired == COMMAND_LINE_PARSER_TRUE) {
                            fprintf(stderr, "%s: Option \"%s\" multiply specified. \n", argv[0], pspec->long_option);
                            return COMMAND_LINE_PARSER_RESULT_OPTION_MULTIPLY_SPECIFIED;
                        }
                        if (pspec->need_argument == COMMAND_LINE_PARSER_TRUE) {
                            /* If the option takes an argument, just go ahead and get the argument */
                            if ((count + 1) == argc) {
                                /* End reached */
                                fprintf(stderr, "%s: Option \"%s\" needs argument. \n", argv[0], pspec->long_option);
                                return COMMAND_LINE_PARSER_RESULT_NOT_SPECIFY_ARGUMENT_TO_OPTION;
                            } else if ((strncmp(argv[count + 1], "--", 2) == 0) || argv[count + 1][0] == '-') {
                                /* Contains other options */
                                /* (Option argument strings beginning with "--" or '-' are not allowed) */
                                fprintf(stderr, "%s: Option \"%s\" needs argument. \n", argv[0], pspec->long_option);
                                return COMMAND_LINE_PARSER_RESULT_NOT_SPECIFY_ARGUMENT_TO_OPTION;
                            }
                            /* Get the option string and move to the next argument string */
                            count++;
                            pspec->argument_string = argv[count];
                        }
                    } else if (arg_str[2 + long_option_len] == '=') {
                        if (pspec->need_argument != COMMAND_LINE_PARSER_TRUE) {
                            /* May be an option containing '='... */
                            continue;
                        }
                        /* An option that has already been acquired was specified */
                        if (pspec->acquired == COMMAND_LINE_PARSER_TRUE) {
                            fprintf(stderr, "%s: Option \"%s\" multiply specified. \n", argv[0], pspec->long_option);
                            return COMMAND_LINE_PARSER_RESULT_OPTION_MULTIPLY_SPECIFIED;
                        }
                        /* Get the option string */
                        pspec->argument_string = &arg_str[2 + long_option_len + 1];
                    } else {
                        /* A longer string is specified. It may match other options, so it will be skipped. */
                        continue;
                    }
                    /* Set to acquired state */
                    pspec->acquired = COMMAND_LINE_PARSER_TRUE;
                    break;
                }
            }
            /* No options found */
            if (spec_no == num_specs) {
                fprintf(stderr, "%s: Unknown long option - \"%s\" \n", argv[0], &arg_str[2]);
                return COMMAND_LINE_PARSER_RESULT_UNKNOWN_OPTION;
            }
        } else if (arg_str[0] == '-') {
            /* Short option (series) */
            uint32_t str_index;
            for (str_index = 1; arg_str[str_index] != '\0'; str_index++) {
                for (spec_no = 0; spec_no < num_specs; spec_no++) {
                    struct CommandLineParserSpecification* pspec = &clps[spec_no];
                    if (arg_str[str_index] == pspec->short_option) {
                        /* An option that has already been acquired was specified */
                        if (pspec->acquired == COMMAND_LINE_PARSER_TRUE) {
                            fprintf(stderr, "%s: Option \'%c\' multiply specified. \n", argv[0], pspec->short_option);
                            return COMMAND_LINE_PARSER_RESULT_OPTION_MULTIPLY_SPECIFIED;
                        }
                        /* Option with arguments */
                        if (pspec->need_argument == COMMAND_LINE_PARSER_TRUE) {
                            /* If the option takes an argument, just go ahead and get the argument */
                            if (arg_str[str_index + 1] != '\0') {
                                /* When taking arguments, the currently selected option must be the last one */
                                fprintf(stderr, "%s: Option \'%c\' needs argument. "
                                        "Please specify tail of short option sequence.\n", argv[0], pspec->short_option);
                                return COMMAND_LINE_PARSER_RESULT_INVAILD_SHORT_OPTION_ARGUMENT;
                            }
                            if ((count + 1) == argc) {
                                /* End reached */
                                fprintf(stderr, "%s: Option \'%c\' needs argument. \n", argv[0], pspec->short_option);
                                return COMMAND_LINE_PARSER_RESULT_NOT_SPECIFY_ARGUMENT_TO_OPTION;
                            } else if ((strncmp(argv[count + 1], "--", 2) == 0) || argv[count + 1][0] == '-') {
                                /* Contains other options */
                                /* ("--" and arguments starting with '-' are not accepted) */
                                fprintf(stderr, "%s: Option \'%c\' needs argument. \n", argv[0], pspec->short_option);
                                return COMMAND_LINE_PARSER_RESULT_NOT_SPECIFY_ARGUMENT_TO_OPTION;
                            }
                            /* Get the option string and move to the next argument string */
                            count++;
                            pspec->argument_string = argv[count];
                        }
                        /* Set to acquired state */
                        pspec->acquired = COMMAND_LINE_PARSER_TRUE;
                        break;
                    }
                }
                /* No options found */
                if (spec_no == num_specs) {
                    fprintf(stderr, "%s: Unknown short option - \'%c\' \n", argv[0], arg_str[str_index]);
                    return COMMAND_LINE_PARSER_RESULT_UNKNOWN_OPTION;
                }
            }
        } else {
            /* A string that is neither an option nor an option argument */
            if (other_string_array == NULL) {
                /* Can't access buffer */
                return COMMAND_LINE_PARSER_RESULT_INVALID_ARGUMENT;
            } else if (other_string_index >= other_string_array_size) {
                /* Buffer size is insufficient */
                fprintf(stderr, "%s: Too many strings specified. \n", argv[0]);
                return COMMAND_LINE_PARSER_RESULT_INSUFFICIENT_OTHER_STRING_ARRAY_SIZE;
            }
            /* Get string */
            other_string_array[other_string_index] = arg_str;
            other_string_index++;
        }
    }

    return COMMAND_LINE_PARSER_RESULT_OK;
}
