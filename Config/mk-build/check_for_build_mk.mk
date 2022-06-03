ifndef BUILD_MK_H__
  ifndef IGNORE_BUILD_MK
    $(info ** You should include build.mk in your arch.make)
    $(info ** to take advantage of the building rules coded in it.)
    $(info ** ---> See top of $(MAIN_OBJDIR)/build.mk )
    $(info **   You can still go ahead using a legacy arch.make,)
    $(info **   as long as you define the proper symbols.)
    $(info **   Use 'make IGNORE_BUILD_MK=1' to bypass this check.)
    $(error --- )
  endif
else
    $(info ** Using build.mk )
endif
