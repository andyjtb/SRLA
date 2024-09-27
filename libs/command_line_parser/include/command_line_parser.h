#ifndef COMMAND_LINE_PARSER_H_INCLDED
#define COMMAND_LINE_PARSER_H_INCLDED

#include <stdint.h>

/* Acquisition result */
typedef enum CommandLineParserResultTag {
    COMMAND_LINE_PARSER_RESULT_OK,                                    /* normal termination */
    COMMAND_LINE_PARSER_RESULT_INVALID_ARGUMENT,                      /* Invalid argument */
    COMMAND_LINE_PARSER_RESULT_INSUFFICIENT_OTHER_STRING_ARRAY_SIZE,  /* The size of the array containing other strings is insufficient */
    COMMAND_LINE_PARSER_RESULT_NOT_SPECIFY_ARGUMENT_TO_OPTION,        /* No argument was specified for an option that requires an argument */
    COMMAND_LINE_PARSER_RESULT_UNKNOWN_OPTION,                        /* An option not defined was specified */
    COMMAND_LINE_PARSER_RESULT_OPTION_MULTIPLY_SPECIFIED,             /* Option specified multiple times */
    COMMAND_LINE_PARSER_RESULT_INVALID_SPECIFICATION,                 /* Invalid specification */
    COMMAND_LINE_PARSER_RESULT_INVAILD_SHORT_OPTION_ARGUMENT          /* Short option argument specified incorrectly */
} CommandLineParserResult;

/* logical constant */
typedef enum CommandLineParserBoolTag {
    COMMAND_LINE_PARSER_FALSE = 0,	/* False */
    COMMAND_LINE_PARSER_TRUE		    /* true */
} CommandLineParserBool;

/* Command line parser specification */
/* Supplementary Note) In the last element of the command line parser specification array, specify 0 for the short option and NULL for the long option. */
struct CommandLineParserSpecification {
    char 				          short_option;		  /* [in] Short option string */
    const char* 		      long_option;		  /* [in] Long option string */
    const char* 		      description;		  /* [in] Argument description */
    CommandLineParserBool	need_argument;		/* [in] Does the option require an argument? */
    const char*				    argument_string;	/* [in,out] Obtained string */
    CommandLineParserBool	acquired;		      /* [out] Was an option specified? */
};

#ifdef __cplusplus
extern "C" {
#endif

/* Print argument description */
void CommandLineParser_PrintDescription(
        const struct CommandLineParserSpecification* clps);

/* Get whether the option was specified from the option name */
CommandLineParserBool CommandLineParser_GetOptionAcquired(
        const struct CommandLineParserSpecification* clps, const char* option_name);

/* Get the option arguments from the option name */
const char* CommandLineParser_GetArgumentString(
        const struct CommandLineParserSpecification* clps, const char* option_name);

/* Parse arguments */
CommandLineParserResult CommandLineParser_ParseArguments(
        struct CommandLineParserSpecification* clps,
        int32_t argc, const char* const* argv,
        const char** other_string_array, uint32_t other_string_array_size);

#ifdef __cplusplus
}
#endif

#endif /* COMMAND_LINE_PARSER_H_INCLDED */
