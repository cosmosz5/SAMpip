import numpy as np
import js_oifits as oifits
import importlib
#importlib.reload(oifits)
import pdb

def oi_merge_wave(filenames):

    ref0 = oifits.open(filenames[0])

    ref0_array = ref0.array
    ref0_flux = ref0.flux
    ref0_wave  = ref0.wavelength
    ref0_v = ref0.vis
    ref0_v2 = ref0.vis2
    ref0_cp = ref0.t3

    eff_wave_f = ref0_wave.eff_wave
    eff_band_f = ref0_wave.eff_band
    
    af_arrayz = ref0_array.arrayz
    af_arrayy = ref0_array.arrayy
    af_arrayx = ref0_array.arrayx
    af_sta_name = ref0_array.sta_name
    af_sta_index = ref0_array.sta_index
    af_diameter = ref0_array.diameter
    af_oirevn = ref0_array.oirevn
    af_frame = ref0_array.frame
    af_tel_name = ref0_array.tel_name
    af_staxyz = ref0_array.staxyz
    af_arrname = ref0_array.arrname

    if ref0_flux != []:
        target_id_f = ref0_flux.target_id
        mjd_f = ref0_flux.mjd
        int_time_f = ref0_flux.int_time
        flux_f = ref0_flux.flux
        fluxerr_f = ref0_flux.fluxerr
        sta_index_f = ref0_flux.sta_index
        flag_f = ref0_flux.flag

    if ref0_v != []:
        target_id_v = ref0_v.target_id
        time_v = ref0_v.time
        mjd_v = ref0_v.mjd
        int_time_v = ref0_v.int_time
        visamp_v = ref0_v.visamp
        visamperr_v = ref0_v.visamperr
        visphi_v = ref0_v.visphi
        visphierr_v = ref0_v.visphierr
        ucoord_v = ref0_v.ucoord
        vcoord_v = ref0_v.vcoord
        sta_index_v = ref0_v.sta_index
        flag_v = ref0_v.flag

    target_id_v2 = ref0_v2.target_id
    time_v2 = ref0_v2.time
    mjd_v2 = ref0_v2.mjd
    int_time_v2 = ref0_v2.int_time
    vis2data_v2 = ref0_v2.vis2data
    vis2err_v2 = ref0_v2.vis2err
    ucoord_v2 = ref0_v2.ucoord
    vcoord_v2 = ref0_v2.vcoord
    sta_index_v2 = ref0_v2.sta_index
    flag_v2 = ref0_v2.flag

    target_id_cp = ref0_cp.target_id
    time_cp = ref0_cp.time
    mjd_cp = ref0_cp.mjd
    int_time_cp = ref0_cp.int_time
    t3amp_cp = ref0_cp.t3amp
    t3amperr_cp = ref0_cp.t3amperr
    t3phi_cp = ref0_cp.t3phi
    t3phierr_cp = ref0_cp.t3phierr
    u1coord_cp = ref0_cp.u1coord
    v1coord_cp = ref0_cp.v1coord
    u2coord_cp = ref0_cp.u2coord
    v2coord_cp = ref0_cp.v2coord
    sta_index_cp = ref0_cp.sta_index
    flag_cp = ref0_cp.flag

    tels_tot = np.array([''.join(ref0_array.sta_name)])

    for i in range(1,len(filenames)):
        print (filenames[i])
        add0 = oifits.open(filenames[i])
        add0_wave = add0.wavelength
        add0_array = add0.array
        add0_flux = add0.flux
        add0_v = add0.vis
        add0_v2 = add0.vis2
        add0_cp = add0.t3

        eff_wave_f = np.append(eff_wave_f, add0_wave.eff_wave)
        eff_band_f = np.append(eff_band_f, add0_wave.eff_band)

        #### Check the different arrays and combine them #####

        tels1 = np.array([''.join(add0_array.sta_name)])
        [index_array] = np.where(tels_tot == tels1)
        if len(index_array) == 0:
            tels_tot = np.append(tels_tot, tels1, axis=0)
            # tels_tot  = tels_tot[0,0,:]
            af_sta_name = np.append(af_sta_name, add0_array.sta_name)
            af_sta_index = np.append(af_sta_index, add0_array.sta_index)
            af_diameter = np.append(af_diameter, add0_array.diameter)
            af_tel_name = np.append(af_tel_name, add0_array.tel_name)
            af_staxyz = np.append(af_staxyz, add0_array.staxyz, axis=0)

        if ref0_flux != []:
            ##### APPEND FLUX #########
            target_id_f = np.append(target_id_f, add0_flux.target_id)
            mjd_f = np.append(mjd_f, add0_flux.mjd)
            int_time_f = np.append(int_time_f, add0_flux.int_time)
            flux_f = np.append(flux_f, add0_flux.flux, axis=0)
            fluxerr_f = np.append(fluxerr_f, add0_flux.fluxerr, axis=0)
            sta_index_f = np.append(sta_index_f, add0_flux.sta_index, axis=0)
            flag_f = np.append(flag_f, add0_flux.flag, axis=0)

        if ref0_v != []:
            ##### APPEND VISIBILITIES #########
            #target_id_v = np.append(target_id_v, )
            #time_v = np.append(time_v, add0_v.time)
            #mjd_v = np.append(mjd_v, add0_v.mjd)
            #int_time_v = np.append(int_time_v, add0_v.int_time)
            visamp_v = np.dstack((visamp_v, add0_v.visamp))
            visamperr_v = np.dstack((visamperr_v, add0_v.visamperr))
            visphi_v = np.dstack((visphi_v, add0_v.visphi))
            visphierr_v = np.dstack((visphierr_v, add0_v.visphierr))
            #ucoord_v = np.dstack((ucoord_v, add0_v.ucoord))
            #vcoord_v = np.dstack((vcoord_v, add0_v.vcoord))
            #sta_index_v = np.dstack((sta_index_v, add0_v.sta_index))
            flag_v = np.dstack((flag_v, add0_v.flag))
            
        ##### APPEND SQUARED VISIBILITIES #########
        #target_id_v2 = np.append(target_id_v2, add0_v2.target_id)
        #time_v2 = np.append(time_v2, add0_v2.time)
        #mjd_v2 = np.append(mjd_v2, add0_v2.mjd)
        #int_time_v2 = np.append(int_time_v2, add0_v2.int_time)
        vis2data_v2 = np.dstack((vis2data_v2, add0_v2.vis2data))
        vis2err_v2 = np.dstack((vis2err_v2, add0_v2.vis2err))
        #ucoord_v2 = np.dstack((ucoord_v2, add0_v2.ucoord))
        #vcoord_v2 = np.dstack((vcoord_v2, add0_v2.vcoord))
        #sta_index_v2 = np.dstack((sta_index_v2, add0_v2.sta_index))
        flag_v2 = np.dstack((flag_v2, add0_v2.flag))
        ##### APPEND CLOSURE PHASES #########
        #target_id_cp = np.append(target_id_cp, add0_cp.target_id)
        #time_cp = np.append(time_cp, add0_cp.time)
        #mjd_cp = np.append(mjd_cp, add0_cp.mjd)
        #int_time_cp = np.append(int_time_cp, add0_cp.int_time)
        t3amp_cp = np.dstack((t3amp_cp, add0_cp.t3amp))
        t3amperr_cp = np.dstack((t3amperr_cp, add0_cp.t3amperr))
        t3phi_cp = np.dstack((t3phi_cp, add0_cp.t3phi))
        t3phierr_cp = np.dstack((t3phierr_cp, add0_cp.t3phierr))
        #u1coord_cp = np.dstack((u1coord_cp, add0_cp.u1coord))
        #v1coord_cp = np.dstack((v1coord_cp, add0_cp.v1coord))
        #u2coord_cp = np.dstack((u2coord_cp, add0_cp.u2coord))
        #v2coord_cp = np.dstack((v2coord_cp, add0_cp.v2coord))
        #sta_index_cp = np.dstack((sta_index_cp, add0_cp.sta_index))
        flag_cp = np.dstack((flag_cp, add0_cp.flag))

    merged = oifits.oifits()
    merged.array = oifits.OI_ARRAY(af_oirevn, af_arrname, af_frame, af_arrayx, af_arrayy, af_arrayz, af_tel_name, \
                                   af_sta_name, af_sta_index, af_diameter, af_staxyz)
    merged.target = ref0.target
    merged.wavelength = oifits.OI_WAVELENGTH(ref0_wave.oirevn, ref0_wave.insname, eff_wave_f, eff_band_f)

    if ref0_flux != []:
        merged.flux = oifits.OI_FLUX(ref0_flux.oirevn, ref0_flux.dateobs, ref0_flux.arrname, ref0_flux.insname,
                                     ref0_flux.calstat, target_id_f, mjd_f, \
                                     int_time_f, flux_f, fluxerr_f, sta_index_f, flag_f)

    if ref0_v != []:
        merged.vis = oifits.OI_VIS(ref0_v.oirevn, ref0_v.dateobs, ref0_v.arrname, ref0_v.insname, target_id_v, time_v, mjd_v, \
                                 int_time_v, np.squeeze(visamp_v), np.squeeze(visamperr_v), np.squeeze(visphi_v), np.squeeze(visphierr_v), \
                                       ucoord_v, vcoord_v, sta_index_v, np.squeeze(flag_v))
    merged.vis2 = oifits.OI_VIS2(ref0_v2.oirevn, ref0_v2.dateobs, ref0_v2.arrname, ref0_v2.insname, target_id_v2, time_v2, mjd_v2, \
                                 int_time_v2, np.squeeze(vis2data_v2), np.squeeze(vis2err_v2), ucoord_v2, vcoord_v2, sta_index_v2, np.squeeze(flag_v2))
    merged.t3 = oifits.OI_T3(ref0_cp.oirevn, ref0_cp.dateobs, ref0_cp.arrname, ref0_cp.insname, target_id_cp, time_cp, mjd_cp, \
                             int_time_cp, np.squeeze(t3amp_cp), np.squeeze(t3amperr_cp), np.squeeze(t3phi_cp), np.squeeze(t3phierr_cp), \
                                 u1coord_cp, v1coord_cp, u2coord_cp, v2coord_cp, sta_index_cp, np.squeeze(flag_cp))

    return merged
