AC_DEFUN([CCS_PPU_EXTRA_FLAGS], [
    AC_ARG_VAR([PPU_EXTRA_CPPFLAGS], [Extra pre-processor flags])
    AC_ARG_VAR([PPU_EXTRA_LDFLAGS], [Extra link flags])
    AC_ARG_VAR([PPU_EXTRA_LIBS], [Extra library flags])

    AC_SUBST(PPU_EXTRA_CPPFLAGS, $PPU_EXTRA_CPPFLAGS)
    AC_SUBST(PPU_EXTRA_LDFLAGS, $PPU_EXTRA_LDFLAGS)
    AC_SUBST(PPU_EXTRA_LIBS, $PPU_EXTRA_LIBS)
])