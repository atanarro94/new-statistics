INC_SRC := ./inc_src
STAT_SRC := ./stats
$(OBJDIR)/statistics.o        :$(STAT_SRC)/statistics.f;                         $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/statistics_2DIO.o      :$(STAT_SRC)/statistics_2DIO.f;                       $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/statistics_2D_usr.o    :$(STAT_SRC)/statistics_2D_usr.f;                     $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/statistics_2D_debug.o    :$(STAT_SRC)/statistics_2D_debug.f;                     $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/time_series.o          :$(STAT_SRC)/time_series.f;                           $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/time_seriesIO.o        :$(STAT_SRC)/time_seriesIO.f;                         $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/pts_redistribute.o     :$(STAT_SRC)/pts_redistribute.f;                      $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/chkpoint.o     :./chkpoint.f;                      $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/IO_tools.o     :./IO_tools.f;                      $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/tripf.o     :./tripf.f;                      $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@


