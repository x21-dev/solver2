/// The header held to the library, when the tests are built: include/solver2.h
/// is read by ImportC (through test/source/solver2_h.c) and compared with
/// what the library has, declaration by declaration, in both directions.
///
/// - Every constant has the same value.
/// - Every struct has the same fields under the same names, of the same
///   types, at the same offsets, and the same size.
/// - Every function and function type has the same return type and the
///   same parameters, names included, and every function is exported.
/// - Neither side has a `solver2_` name the other lacks.
///
/// A mismatch fails the build of the tests and names the declaration. Only
/// the tests need a C compiler (ImportC runs its preprocessor); the library
/// builds without one.
module solver2c.conformance;

version (unittest):

static import api = solver2c.api;
static import d = solver2c.c;
static import h = solver2_h;

/// Whether `name` is one of the interface's.
private bool ours(string name) pure nothrow @nogc @safe
{
    return name.length > 8 && (name[0 .. 8] == "solver2_" || name[0 .. 8] == "SOLVER2_");
}

/// A type's text as both sides can be compared by: without `const`, which
/// ImportC keeps on a pointer to a number and drops on a pointer to a
/// struct, and without the attributes D's own functions carry. Neither
/// changes how a call is made.
private string plain(string type) pure @safe
{
    string o;
    bool[] wrappers;                // per open bracket: is it a const( ... )?
    static immutable string[] attributes = [" nothrow", " @nogc", " @system", " @trusted", " @safe", " pure"];
    next: for (size_t i = 0; i < type.length;)
    {
        foreach (word; attributes)
            if (type[i .. $].length >= word.length && type[i .. i + word.length] == word)
            {
                i += word.length;
                continue next;
            }
        if (type[i .. $].length >= 6 && type[i .. i + 6] == "const(")
        {
            wrappers ~= true;
            i += 6;
            continue;
        }
        if (type[i] == '(')
            wrappers ~= false;
        else if (type[i] == ')')
        {
            immutable bool wrapper = wrappers[$ - 1];
            wrappers = wrappers[0 .. $ - 1];
            if (wrapper)
            {
                ++i;
                continue;
            }
        }
        o ~= type[i++];
    }
    return o;
}

static assert(plain("extern (C) int function(const(solver2_gene)* genes, const(char*)* names) nothrow @nogc")
              == "extern (C) int function(solver2_gene* genes, char** names)");

// The header's declarations, each against the library's.
static foreach (name; __traits(allMembers, h))
    static if (ours(name))
    {
        static if (is(__traits(getMember, h, name) == struct))
        {
            static assert(is(__traits(getMember, d, name) == struct), "solver2.h's struct " ~ name ~ " is not in c.d");
            static if (__traits(getMember, h, name).tupleof.length > 0)
            {
                static assert(__traits(getMember, h, name).tupleof.length == __traits(getMember, d, name).tupleof.length
                              && __traits(getMember, h, name).sizeof == __traits(getMember, d, name).sizeof,
                              "struct " ~ name ~ ": another number of fields, or another size, in solver2.h than in c.d");
                static foreach (i; 0 .. __traits(getMember, h, name).tupleof.length)
                    static assert(__traits(identifier, __traits(getMember, h, name).tupleof[i])
                                  == __traits(identifier, __traits(getMember, d, name).tupleof[i])
                                  && __traits(getMember, h, name).tupleof[i].offsetof
                                  == __traits(getMember, d, name).tupleof[i].offsetof
                                  && plain(typeof(__traits(getMember, h, name).tupleof[i]).stringof)
                                  == plain(typeof(__traits(getMember, d, name).tupleof[i]).stringof),
                                  "struct " ~ name ~ ": field " ~ __traits(identifier, __traits(getMember, h, name).tupleof[i])
                                  ~ " of solver2.h is another in c.d (its name, its type or its place)");
            }
        }
        else static if (is(typeof(__traits(getMember, h, name)) == function))
        {
            static assert(__traits(hasMember, api, name), "solver2.h's function " ~ name ~ " is not in api.d");
            static assert(__traits(getVisibility, __traits(getMember, api, name)) == "export",
                          name ~ " is not exported from the library");
            static assert(plain(typeof(&__traits(getMember, h, name)).stringof)
                          == plain(typeof(&__traits(getMember, api, name)).stringof),
                          name ~ ": solver2.h declares " ~ plain(typeof(&__traits(getMember, h, name)).stringof)
                          ~ ", api.d defines " ~ plain(typeof(&__traits(getMember, api, name)).stringof));
        }
        else static if (is(__traits(getMember, h, name)))
        {
            static assert(__traits(hasMember, d, name), "solver2.h's type " ~ name ~ " is not in c.d");
            static assert(plain(__traits(getMember, h, name).stringof) == plain(__traits(getMember, d, name).stringof),
                          name ~ ": solver2.h has " ~ plain(__traits(getMember, h, name).stringof) ~ ", c.d has "
                          ~ plain(__traits(getMember, d, name).stringof));
        }
        else
        {
            static assert(__traits(hasMember, d, name), "solver2.h's constant " ~ name ~ " is not in c.d");
            static assert(cast(long) __traits(getMember, h, name) == cast(long) __traits(getMember, d, name),
                          name ~ " has another value in solver2.h than in c.d");
        }
    }

// ... and nothing of the library's that the header lacks.
static foreach (name; __traits(allMembers, d))
    static if (ours(name))
        static assert(__traits(hasMember, h, name), "c.d's " ~ name ~ " is not in solver2.h");

static foreach (name; __traits(allMembers, api))
    static if (ours(name))
        static assert(__traits(hasMember, h, name), "api.d's " ~ name ~ " is not in solver2.h");

unittest
{
    // The comparison itself has something to compare: the header came
    // through, with its functions and its structs.
    size_t functions, structs;
    static foreach (name; __traits(allMembers, h))
        static if (ours(name))
        {
            static if (is(typeof(__traits(getMember, h, name)) == function))
                ++functions;
            else static if (is(__traits(getMember, h, name) == struct))
                ++structs;
        }
    assert(functions == 37 && structs == 18);
}
