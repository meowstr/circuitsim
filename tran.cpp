#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <unistd.h>

#include <vector>

const char * delims = " \n\t";

struct capacitor_t {
    char * name;
    char * n1;
    char * n2;
    double farads;
    double accum_volts;
    int line;
};

struct inductor_t {
    char * name;
    char * n1;
    char * n2;
    double henries;
    double accum_amps;
    int line;
};

static std::vector< char * > line_list;
static std::vector< capacitor_t > capacitor_list;
static std::vector< inductor_t > inductor_list;

static std::vector< char * > output_name_list;
static std::vector< double > output_value_list;

static void parse_c( char * command )
{
    char * n1 = strtok( nullptr, delims );
    char * n2 = strtok( nullptr, delims );
    char * farads_str = strtok( nullptr, delims );

    capacitor_t c;
    c.name = strdup( command );
    c.n1 = strdup( n1 );
    c.n2 = strdup( n2 );
    c.farads = atof( farads_str );
    c.accum_volts = 0;
    c.line = line_list.size();

    capacitor_list.push_back( c );
}

static void parse_l( char * command )
{
    char * n1 = strtok( nullptr, delims );
    char * n2 = strtok( nullptr, delims );
    char * henries_str = strtok( nullptr, delims );

    inductor_t l;
    l.name = strdup( command );
    l.n1 = strdup( n1 );
    l.n2 = strdup( n2 );
    l.henries = atof( henries_str );
    l.accum_amps = 0;
    l.line = line_list.size();

    inductor_list.push_back( l );
}

static double find_output_param( char * name )
{
    for ( int i = 0; i < output_name_list.size(); i++ ) {
        if ( strcmp( output_name_list[ i ], name ) == 0 ) {
            return output_value_list[ i ];
        }
    }

    printf( "missing output param: %s\n", name );
    return 0;
}

static void update_capacitor( capacitor_t & c, double dt )
{
    char buffer[ 1024 ];
    snprintf( buffer, 1024, "i_v_%s", c.name );
    double i = find_output_param( buffer );

    c.accum_volts += dt * -i / c.farads;
}

static void update_inductor( inductor_t & l, double dt )
{
    double v1 = find_output_param( l.n1 );
    double v2 = find_output_param( l.n2 );

    // maybe wrong, but i hate inductors so eh
    l.accum_amps += dt * ( v1 - v2 ) / l.henries;
}

static void update_components( double dt )
{
    for ( capacitor_t & c : capacitor_list ) {
        update_capacitor( c, dt );
    }

    for ( inductor_t & l : inductor_list ) {
        update_inductor( l, dt );
    }
}

static void iterate( double dt )
{
    for ( capacitor_t & c : capacitor_list ) {
        snprintf(
            line_list[ c.line ],
            1024,
            "v_%s %s %s %le",
            c.name,
            c.n1,
            c.n2,
            c.accum_volts
        );
    }

    for ( inductor_t & l : inductor_list ) {
        snprintf(
            line_list[ l.line ],
            1024,
            "i_%s %s %s %le",
            l.name,
            l.n1,
            l.n2,
            l.accum_amps
        );
    }

    int parent_to_child[ 2 ];
    int child_to_parent[ 2 ];
    pipe( parent_to_child );
    pipe( child_to_parent );

    pid_t pid = fork();

    if ( pid == 0 ) {
        dup2( parent_to_child[ 0 ], STDIN_FILENO );
        dup2( child_to_parent[ 1 ], STDOUT_FILENO );
        close( parent_to_child[ 0 ] );
        close( parent_to_child[ 1 ] );
        close( child_to_parent[ 0 ] );
        close( child_to_parent[ 1 ] );

        execl( "./sim", "sim", "-s", nullptr );

        printf( "spawning sim process failed\n" );

    } else {
        close( parent_to_child[ 0 ] );
        close( child_to_parent[ 1 ] );

        for ( char * line : line_list ) {
            write( parent_to_child[ 1 ], line, strlen( line ) );
            write( parent_to_child[ 1 ], "\n", 1 );
        }

        close( parent_to_child[ 1 ] );

        FILE * f = fdopen( child_to_parent[ 0 ], "r" );

        // output_name_list.clear();
        output_value_list.clear();

        int fill_names = output_name_list.empty();

        char line[ 1024 ];
        while ( fgets( line, 1024, f ) ) {
            char * name = strtok( line, delims );
            char * equal_sign = strtok( nullptr, delims );
            char * value_str = strtok( nullptr, delims );
            if ( fill_names ) output_name_list.push_back( strdup( name ) );
            output_value_list.push_back( atof( value_str ) );
        }
        fclose( f );

        update_components( dt );
    }
}

int main( int argc, char ** argv )
{
    double dt = 0.1;
    double steps = 100;
    int plot_enabled = 0;

    std::vector< char * > plot_col_names;

    FILE * f = stdin;

    char line[ 1024 ];
    char buffer[ 1024 ];
    while ( fgets( line, 1024, f ) ) {
        strncpy( buffer, line, 1024 );
        char * command = strtok( buffer, delims );

        if ( !command ) continue;

        char * new_line = (char *) malloc( 1024 );
        if ( command[ 0 ] == 'c' )
            parse_c( command );
        else if ( command[ 0 ] == 'l' )
            parse_l( command );
        else if ( strcmp( command, ".step" ) == 0 ) {
            char * step_size_str = strtok( nullptr, delims );
            char * step_count_str = strtok( nullptr, delims );
            dt = atof( step_size_str );
            steps = atof( step_count_str );
        } else if ( strcmp( command, ".plot" ) == 0 ) {
            for ( ;; ) {
                char * name = strtok( nullptr, delims );
                if ( name )
                    plot_col_names.push_back( strdup( name ) );
                else
                    break;
            }
            plot_enabled = 1;
        } else
            strncpy( new_line, line, 1024 );

        line_list.push_back( new_line );
    }

    double time = 0;
    iterate( dt ); // first iteration to get the column names

    FILE * out;
    if ( plot_enabled ) {
        out = popen( "gnuplot -persist", "w" );
        if ( !out ) {
            printf( "failed to open gnuplot\n" );
            return 1;
        }

        int col = 1;
        for ( int i = 0; i < output_name_list.size(); i++ ) {
            if ( strcmp( output_name_list[ i ], plot_col_names[ 0 ] ) == 0 )
                col = i + 2;
        }

        fprintf( out, "set xlabel 'time (s)'\n" );
        fprintf(
            out,
            "set print $DATA; lines = system('cat'); print lines; unset print; "
            "plot $DATA using 1:%d with lines title 'node \"%s\"'\n",
            col,
            plot_col_names[ 0 ]
        );
    } else {
        out = stdout;
    }

    for ( int i = 0; i < steps; i++ ) {
        if ( !plot_enabled && i == 0 ) {
            fprintf( out, "time " );
            for ( char * name : output_name_list ) {
                fprintf( out, "%s ", name );
            }
            fprintf( out, "\n" );
        }

        fprintf( out, "% le ", time );
        for ( double value : output_value_list ) {
            fprintf( out, "% le ", value );
        }
        fprintf( out, "\n" );

        iterate( dt );
        time += dt;
    }

    fflush( out );
    pclose( out );

    return 0;
}
