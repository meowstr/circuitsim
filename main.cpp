// g++ main.cpp -llapacke

#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> // getopt

#include <vector>

static const char * delims = " \n\t";

static void evaluate_system( double * x, double * out );
static void evaluate_jacobian( double * x, double * out );

static int verbose = 0;

static void root_newton_step( int n, double * x )
{
    static int old_n = 0;

    static double * a;
    static double * b;
    static int * ipiv;

    static double * fx;
    static double * old_fx;
    static double * new_x;

    if ( n != old_n ) {
        if ( verbose ) printf( "root_newton_step: allocate (n = %d)\n", n );
        a = (double *) malloc( n * n * sizeof( double ) );
        b = (double *) malloc( n * sizeof( double ) );
        ipiv = (int *) malloc( n * sizeof( int ) );

        fx = (double *) malloc( n * sizeof( double ) );
        old_fx = (double *) malloc( n * sizeof( double ) );
        new_x = (double *) malloc( ( n + 1 ) * sizeof( double ) );

        old_n = n;
    }

    // setup A * X = B  =>  J(x) * (new_x - x) = -F(x)

    evaluate_jacobian( x, a );
    evaluate_system( x, b );

    for ( int i = 0; i < n; i++ )
        old_fx[ i ] = b[ i ];
    for ( int i = 0; i < n; i++ )
        b[ i ] = -b[ i ];

    int info = LAPACKE_dgesv( LAPACK_ROW_MAJOR, n, 1, a, n, ipiv, b, 1 );

    // solution stored in b

    if ( info ) {
        if ( verbose ) printf( "root_newton_step: dgesv failed\n" );
        return;
    }

    double alpha = 1;
backtrack:
    new_x[ 0 ] = x[ 0 ];
    for ( int i = 0; i < n; i++ ) {
        new_x[ i + 1 ] = x[ i + 1 ] + alpha * b[ i ];
    }
    evaluate_system( new_x, fx );

    // keep changes in x small
    double norm1 = 0;
    for ( int i = 0; i < n + 1; i++ )
        norm1 = fmax( norm1, abs( new_x[ i ] - x[ i ] ) );
    // if ( norm1 > 1. ) {
    //     alpha /= 2.;
    //     goto backtrack;
    // }

    // keep changes in fx small
    double norm2 = 0;
    for ( int i = 0; i < n; i++ )
        norm2 = fmax( norm2, abs( fx[ i ] - old_fx[ i ] ) );
    if ( norm2 > 10. ) {
        alpha /= 2.;
        if ( verbose ) {
            printf(
                "root_newton_step: backtrack (dx = %.2le df = %.2le) \n",
                norm1,
                norm2
            );
        }
        goto backtrack;
    }

    for ( int i = 0; i < n; i++ )
        x[ i + 1 ] = new_x[ i + 1 ];
}

static void root_newton( int n, double * x )
{
    static int old_n = 0;

    static double * fx;

    if ( n != old_n ) {
        fx = (double *) malloc( n * sizeof( double ) );
        old_n = n;
    }

    for ( int i = 0; i < 20; i++ ) {
        if ( verbose ) {
            evaluate_system( x, fx );

            printf( "root_newton: " );
            printf( "F( " );
            for ( int k = 0; k < n; k++ ) {
                printf( "% .2le ", x[ k + 1 ] );
            }
            printf( ") = " );
            for ( int k = 0; k < n; k++ ) {
                printf( "% .2le ", fx[ k ] );
            }
            printf( "\n" );
        }

        root_newton_step( n, x );
    }

    if ( verbose ) printf( "\n" );
}

static double V_A = 10.; // early voltage

static double eber_moll( double vbe, double vce )
{
    return 1E-13 * ( exp( vbe * 40. ) - 1. ) * ( 1. + vce / V_A );
}

static double eber_moll_partial_vbe( double vbe, double vce )
{
    return 40. * eber_moll( vbe, vce );
}

static double eber_moll_partial_vce( double vbe, double vce )
{
    return 1E-13 * ( exp( vbe * 40. ) - 1. ) * ( 1. / V_A );
}

static double early_beta( double vce )
{
    double beta = 100;
    return beta * ( 1 + vce / V_A );
}

static double early_beta_partial()
{
    double beta = 100;
    return beta / V_A;
}

struct equation_t {
    char * name;
    std::vector< char * > term_list;
};

static struct system_t {
    int next_current_id;
    int next_voltage_id;

    std::vector< char * > unknown_list;

    std::vector< equation_t > equation_list;
} sys;

static int find_unknown_index( char * unknown )
{
    for ( int i = 0; i < sys.unknown_list.size(); i++ ) {
        if ( strcmp( sys.unknown_list[ i ], unknown ) == 0 ) return i;
    }
    printf( "failed to find unknown: %s\n", unknown );
    return 0;
}

static double
evaluate_term_partial( double * x, char * term, char * unknown, char * eq_name )
{
    char buffer[ 1024 ];
    strncpy( buffer, term, 1024 );

    char * command = strtok( buffer, delims );

    if ( term[ 0 ] == 'i' ) {
        // if not w.r.t to the current then its 0
        if ( strcmp( command, unknown ) != 0 ) return 0;

        char * sign = strtok( nullptr, delims );
        if ( sign[ 0 ] == '+' ) {
            return 1;
        } else if ( sign[ 0 ] == '-' ) {
            return -1;
        } else {
            printf( "i term missing sign: %s\n", term );
            return 0;
        }
    }

    if ( term[ 0 ] == 'r' ) {
        char * n1_str = strtok( nullptr, delims );
        char * n2_str = strtok( nullptr, delims );
        char * ohms_str = strtok( nullptr, delims );

        int n1 = find_unknown_index( n1_str );
        int n2 = find_unknown_index( n2_str );
        double ohms = atof( ohms_str );

        double di;
        if ( strcmp( n1_str, unknown ) == 0 )
            di = 1. / ohms;
        else if ( strcmp( n2_str, unknown ) == 0 )
            di = -1. / ohms;
        else
            return 0;

        if ( strcmp( n1_str, eq_name ) == 0 ) return -di;
        if ( strcmp( n2_str, eq_name ) == 0 ) return di;

        printf( "resistor not in kcl equation\n" );
        return 0;
    }

    if ( term[ 0 ] == 'd' ) {
        char * n1_str = strtok( nullptr, delims );
        char * n2_str = strtok( nullptr, delims );

        int n1 = find_unknown_index( n1_str );
        int n2 = find_unknown_index( n2_str );

        double di;
        if ( strcmp( n1_str, unknown ) == 0 )
            di = eber_moll_partial_vbe( x[ n1 ] - x[ n2 ], 0 );
        else if ( strcmp( n2_str, unknown ) == 0 )
            di = -eber_moll_partial_vbe( x[ n1 ] - x[ n2 ], 0 );
        else
            return 0;

        if ( strcmp( n1_str, eq_name ) == 0 ) return -di;
        if ( strcmp( n2_str, eq_name ) == 0 ) return di;

        printf( "diode not in kcl equation\n" );
        return 0;
    }

    if ( term[ 0 ] == 'v' ) {
        char * n1_str = strtok( nullptr, delims );
        char * n2_str = strtok( nullptr, delims );
        char * voltage_str = strtok( nullptr, delims );

        int n1 = find_unknown_index( n1_str );
        int n2 = find_unknown_index( n2_str );
        double voltage = atof( voltage_str );

        if ( strcmp( n1_str, unknown ) == 0 ) return -1;
        if ( strcmp( n2_str, unknown ) == 0 ) return 1;
        return 0;
    }

    if ( term[ 0 ] == 'q' ) {
        char * nc_str = strtok( nullptr, delims );
        char * nb_str = strtok( nullptr, delims );
        char * ne_str = strtok( nullptr, delims );
        char * bjt_type = strtok( nullptr, delims );

        int nc = find_unknown_index( nc_str );
        int nb = find_unknown_index( nb_str );
        int ne = find_unknown_index( ne_str );

        // transport equations:
        // https://en.wikipedia.org/wiki/Bipolar_junction_transistor#Ebers%E2%80%93Moll_model
        double Is = 1E-13;
        double Vt = 25. / 1000.;
        double beta_r = 10.;
        double beta_f = 100.;
        double Vbe = x[ nb ] - x[ ne ];
        double Vbc = x[ nb ] - x[ nc ];
        double dir = 1.;

        if ( strcmp( bjt_type, "npn" ) == 0 ) {
            dir = 1.;
        } else if ( strcmp( bjt_type, "pnp" ) == 0 ) {
            dir = 1.; // derivatives flip sign again
            Vbe *= -1.;
            Vbc *= -1.;
        } else {
            printf( "unknown bjt type: %s\n", bjt_type );
            return 0;
        }
        if ( strcmp( eq_name, nc_str ) == 0 ) {
            if ( strcmp( unknown, nc_str ) == 0 )
                return dir * ( -Is / Vt ) *
                       ( exp( Vbc / Vt ) + ( 1. / beta_r ) * exp( Vbc / Vt ) );
            if ( strcmp( unknown, nb_str ) == 0 )
                return dir * ( -Is / Vt ) *
                       ( ( exp( Vbe / Vt ) - exp( Vbc / Vt ) ) -
                         ( 1. / beta_r ) * ( exp( Vbc / Vt ) ) );
            if ( strcmp( unknown, ne_str ) == 0 )
                return dir * -( -Is / Vt ) * exp( Vbe / Vt );

            return 0;

            // return -Is * ( ( exp( Vbe / Vt ) - exp( Vbc / Vt ) ) -
            //               ( 1. / beta_r ) * ( exp( Vbc / Vt ) - 1. ) );
        }
        if ( strcmp( eq_name, nb_str ) == 0 ) {
            if ( strcmp( unknown, nc_str ) == 0 )
                return dir * -( -Is / Vt ) * ( 1. / beta_r ) * exp( Vbc / Vt );
            if ( strcmp( unknown, nb_str ) == 0 )
                return dir * ( -Is / Vt ) *
                       ( ( 1. / beta_f ) * ( exp( Vbe / Vt ) ) +
                         ( 1. / beta_r ) * ( exp( Vbc / Vt ) ) );
            if ( strcmp( unknown, ne_str ) == 0 )
                return dir * -( -Is / Vt ) * ( 1. / beta_f ) * exp( Vbe / Vt );

            return 0;

            // return -Is * ( ( 1. / beta_f ) * ( exp( Vbe / Vt ) - 1. ) +
            //               ( 1. / beta_r ) * ( exp( Vbc / Vt ) - 1. ) );
        }
        if ( strcmp( eq_name, ne_str ) == 0 ) {
            if ( strcmp( unknown, nc_str ) == 0 )
                return dir * ( Is / Vt ) * exp( Vbc / Vt );
            if ( strcmp( unknown, nb_str ) == 0 )
                return dir * ( Is / Vt ) *
                       ( ( exp( Vbe / Vt ) - exp( Vbc / Vt ) ) +
                         ( 1. / beta_f ) * ( exp( Vbe / Vt ) ) );
            if ( strcmp( unknown, ne_str ) == 0 )
                return dir * ( Is / Vt ) *
                       ( -exp( Vbe / Vt ) +
                         ( 1. / beta_f ) * ( -exp( Vbe / Vt ) ) );

            return 0;

            // return Is * ( ( exp( Vbe / Vt ) - exp( Vbc / Vt ) ) +
            //               ( 1. / beta_f ) * ( exp( Vbe / Vt ) - 1. ) );
        }

        printf( "bjt not in kcl equation\n" );
        return 0;
    }

    printf( "unknown equation term: %s\n", term );

    return 0;
}

static double evaluate_partial( double * x, equation_t & e, char * unknown )
{
    double sum = 0;
    for ( char * term : e.term_list ) {
        sum += evaluate_term_partial( x, term, unknown, e.name );
    }
    return sum;
}

static void evaluate_jacobian( double * x, double * out )
{
    int i = 0;
    for ( equation_t & e : sys.equation_list ) {
        for ( char * unknown : sys.unknown_list ) {
            if ( strcmp( unknown, "0" ) == 0 ) continue; // skip the dummy one
            out[ i++ ] = evaluate_partial( x, e, unknown );
        }
    }
}

static double evaluate_term( double * x, char * term, char * eq_name )
{
    char buffer[ 1024 ];
    strncpy( buffer, term, 1024 );

    char * command = strtok( buffer, delims );

    if ( term[ 0 ] == 'i' ) {
        int i = find_unknown_index( command );
        char * sign = strtok( nullptr, delims );
        if ( sign[ 0 ] == '+' ) {
            return x[ i ];
        } else if ( sign[ 0 ] == '-' ) {
            return -x[ i ];
        } else {
            printf( "i term missing sign: %s\n", term );
            return 0;
        }
    }

    if ( term[ 0 ] == 'r' ) {
        char * n1_str = strtok( nullptr, delims );
        char * n2_str = strtok( nullptr, delims );
        char * ohms_str = strtok( nullptr, delims );

        int n1 = find_unknown_index( n1_str );
        int n2 = find_unknown_index( n2_str );
        double ohms = atof( ohms_str );

        double i = ( x[ n1 ] - x[ n2 ] ) / ohms;

        if ( strcmp( n1_str, eq_name ) == 0 ) return -i;
        if ( strcmp( n2_str, eq_name ) == 0 ) return i;

        printf( "resistor not in kcl equation\n" );
        return 0;
    }

    if ( term[ 0 ] == 'd' ) {
        char * n1_str = strtok( nullptr, delims );
        char * n2_str = strtok( nullptr, delims );

        int n1 = find_unknown_index( n1_str );
        int n2 = find_unknown_index( n2_str );

        double i = eber_moll( x[ n1 ] - x[ n2 ], 0 );

        if ( strcmp( n1_str, eq_name ) == 0 ) return -i;
        if ( strcmp( n2_str, eq_name ) == 0 ) return i;

        printf( "diode not in kcl equation\n" );
        return 0;
    }

    if ( term[ 0 ] == 'v' ) {
        char * n1_str = strtok( nullptr, delims );
        char * n2_str = strtok( nullptr, delims );
        char * voltage_str = strtok( nullptr, delims );

        int n1 = find_unknown_index( n1_str );
        int n2 = find_unknown_index( n2_str );
        double voltage = atof( voltage_str );

        return ( x[ n2 ] - x[ n1 ] ) - voltage;
    }

    if ( term[ 0 ] == 'q' ) {
        char * nc_str = strtok( nullptr, delims );
        char * nb_str = strtok( nullptr, delims );
        char * ne_str = strtok( nullptr, delims );
        char * bjt_type = strtok( nullptr, delims );

        int nc = find_unknown_index( nc_str );
        int nb = find_unknown_index( nb_str );
        int ne = find_unknown_index( ne_str );

        // transport equations:
        // https://en.wikipedia.org/wiki/Bipolar_junction_transistor#Ebers%E2%80%93Moll_model
        double Is = 1E-13;
        double Vt = 25. / 1000.;
        double beta_r = 10.;
        double beta_f = 100.;
        double Vbe = x[ nb ] - x[ ne ];
        double Vbc = x[ nb ] - x[ nc ];
        double dir = 1.;

        if ( strcmp( bjt_type, "npn" ) == 0 ) {
            dir = 1.;
        } else if ( strcmp( bjt_type, "pnp" ) == 0 ) {
            dir = -1.;
            Vbe *= -1;
            Vbc *= -1;
        } else {
            printf( "unknown bjt type: %s\n", bjt_type );
            return 0;
        }

        if ( strcmp( eq_name, nc_str ) == 0 ) {
            // NOTE: negative because I_c sinks current
            return dir * -Is *
                   ( ( exp( Vbe / Vt ) - exp( Vbc / Vt ) ) -
                     ( 1. / beta_r ) * ( exp( Vbc / Vt ) - 1. ) );
        }
        if ( strcmp( eq_name, nb_str ) == 0 ) {
            // NOTE: negative because I_b sinks current
            return dir * -Is *
                   ( ( 1. / beta_f ) * ( exp( Vbe / Vt ) - 1. ) +
                     ( 1. / beta_r ) * ( exp( Vbc / Vt ) - 1. ) );
        }
        if ( strcmp( eq_name, ne_str ) == 0 ) {
            return dir * Is *
                   ( ( exp( Vbe / Vt ) - exp( Vbc / Vt ) ) +
                     ( 1. / beta_f ) * ( exp( Vbe / Vt ) - 1. ) );
        }

        printf( "bjt not in kcl equation\n" );
        return 0;
    }

    printf( "unknown equation term: %s\n", term );

    return 0;
}

static double evaluate_equation( double * x, equation_t & e )
{
    double sum = 0;
    for ( char * term : e.term_list ) {
        sum += evaluate_term( x, term, e.name );
    }
    return sum;
}

static void evaluate_system( double * x, double * out )
{
    for ( int i = 0; i < sys.equation_list.size(); i++ ) {
        out[ i ] = evaluate_equation( x, sys.equation_list[ i ] );
    }
}

static void add_to_kcl( char * node, char * term )
{
    // ignore KCL equation for ground
    if ( strcmp( node, "0" ) == 0 ) return;

    equation_t * eq = nullptr;
    for ( equation_t & e : sys.equation_list ) {
        if ( strcmp( e.name, node ) == 0 ) eq = &e;
    }

    if ( eq == nullptr ) {
        equation_t e;
        e.name = strdup( node );
        sys.equation_list.push_back( e );
        sys.unknown_list.push_back( strdup( node ) );
        eq = &sys.equation_list.back();
    }

    eq->term_list.push_back( strdup( term ) );
}

static char buffer[ 1024 ];
static void parse_v( char * command )
{
    char * n1 = strtok( nullptr, delims );
    char * n2 = strtok( nullptr, delims );
    char * voltage = strtok( nullptr, delims );

    // add new current unknown: i
    // add new equation: v2 - v1 = voltage
    // add to n1 KCL: -i
    // add to n2 KCL: i

    char current_name[ 32 ];
    snprintf( current_name, 32, "i_%s", command );
    sys.unknown_list.push_back( strdup( current_name ) );

    snprintf( buffer, 1024, "v %s %s %s", n1, n2, voltage );
    equation_t e;
    e.name = strdup( "vsource" );
    e.term_list.push_back( strdup( buffer ) );
    sys.equation_list.push_back( e );

    snprintf( buffer, 1024, "%s -", current_name );
    add_to_kcl( n1, buffer );

    snprintf( buffer, 1024, "%s +", current_name );
    add_to_kcl( n2, buffer );
}

static void parse_d( char * command )
{
    char * n1 = strtok( nullptr, delims );
    char * n2 = strtok( nullptr, delims );

    // add to n1 KCL: -Is e^((v1 - v2) / V_T)
    // add to n2 KCL: Is e^((v1 - v2) / V_T)

    snprintf( buffer, 1024, "d %s %s", n1, n2 );
    add_to_kcl( n1, buffer );
    add_to_kcl( n2, buffer );
}

static void parse_r( char * command )
{
    char * n1 = strtok( nullptr, delims );
    char * n2 = strtok( nullptr, delims );
    char * ohms = strtok( nullptr, delims );

    // add to n1 KCL: -(v1 - v2) / ohms
    // add to n2 KCL: (v1 - v2) / ohms

    snprintf( buffer, 1024, "r %s %s %s", n1, n2, ohms );
    add_to_kcl( n1, buffer );
    add_to_kcl( n2, buffer );
}

static void parse_q( char * command )
{
    char * nc = strtok( nullptr, delims );
    char * nb = strtok( nullptr, delims );
    char * ne = strtok( nullptr, delims );
    char * type = strtok( nullptr, delims );

    // add transport equations to each kcl

    snprintf( buffer, 1024, "q %s %s %s %s", nc, nb, ne, type );
    add_to_kcl( nc, buffer );
    add_to_kcl( nb, buffer );
    add_to_kcl( ne, buffer );
}

int main( int argc, char ** argv )
{
    verbose = 1;

    int opt;
    while ( ( opt = getopt( argc, argv, "s" ) ) != -1 ) {
        if ( opt == 's' ) {
            verbose = 0;
        }
        if ( opt == '?' ) {
            printf( "%s [-s] filename\n", argv[ 0 ] );
            return 1;
        }
    }

    // dummy unknown for ground, not actually solved for
    sys.unknown_list.push_back( strdup( "0" ) );

    FILE * f;
    if ( optind < argc ) {
        f = fopen( argv[ optind ], "r" );
    } else {
        f = stdin;
        printf( "reading stdin\n" );
    }

    if ( !f ) {
        printf( "failed to read file\n" );
        return 1;
    }

    char line[ 1024 ];
    while ( fgets( line, 1024, f ) ) {
        char * command = strtok( line, delims );

        if ( !command ) continue;

        if ( command[ 0 ] == 'v' )
            parse_v( command );
        else if ( command[ 0 ] == 'r' )
            parse_r( command );
        else if ( command[ 0 ] == 'd' )
            parse_d( command );
        else if ( command[ 0 ] == 'q' )
            parse_q( command );
        else if ( strcmp( command, ".end" ) == 0 )
            break;
        else
            printf( "wtf did u type: %s\n", command );
    }

    fclose( f );

    if ( verbose ) {
        printf( "unknowns: " );
        for ( int i = 1; i < sys.unknown_list.size(); i++ ) {
            printf( "%s ", sys.unknown_list[ i ] );
        }
        printf( "\n" );

        printf( "equations:\n" );
        for ( equation_t & eq : sys.equation_list ) {
            printf( "\tname: %s\n", eq.name );
            for ( char * term : eq.term_list ) {
                printf( "\t\t%s\n", term );
            }
        }
        printf( "\n" );
    }

    int n = sys.equation_list.size();
    double * x = (double *) malloc( ( n + 1 ) * sizeof( double ) );
    double * fx = (double *) malloc( n * sizeof( double ) );
    double * j = (double *) malloc( n * n * sizeof( double ) );

    for ( int i = 0; i < n + 1; i++ )
        x[ i ] = 0;

    // x[ 2 ] = 10;

    evaluate_system( x, fx );

    if ( verbose ) {
        printf( "system at 0: " );
        for ( int i = 0; i < n; i++ )
            printf( "% .2le ", fx[ i ] );
        printf( "\n" );

        evaluate_jacobian( x, j );
        printf( "jacobian at 0: " );
        for ( int i = 0; i < n * n; i++ ) {
            if ( i % n == 0 ) printf( "\n\t" );
            printf( "% .2le ", j[ i ] );
        }
        printf( "\n" );
    }

    root_newton( n, x );

    if ( verbose ) printf( "solve:\n" );
    for ( int i = 0; i < n; i++ ) {
        if ( verbose ) printf( "\t" );
        printf( "%5s = % le\n", sys.unknown_list[ i + 1 ], x[ i + 1 ] );
    }

    return 0;
}
