add_rules("mode.debug", "mode.release")

set_languages("c++17")

target("Ephemeris")
    set_kind("binary")

    add_includedirs(
        "src",
        "src/integrators",
        "src/math",
        "src/physics",
        "src/tests",
        "src/utils",
        "src/modules"
    )

    add_files("src/**.cpp", "src/**.cu")
    add_files("/*.rc")
target_end()