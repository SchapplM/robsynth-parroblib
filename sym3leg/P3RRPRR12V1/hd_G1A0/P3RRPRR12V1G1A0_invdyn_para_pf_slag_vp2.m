% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR12V1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:02
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:01:39
% EndTime: 2020-08-06 19:01:46
% DurationCPUTime: 6.82s
% Computational Cost: add. (22134->456), mult. (35376->760), div. (4953->6), fcn. (32256->18), ass. (0->311)
t221 = pkin(1) + pkin(2);
t205 = sin(qJ(2,3));
t218 = xDP(3);
t219 = xDP(2);
t220 = xDP(1);
t223 = 0.1e1 / qJ(3,3);
t211 = cos(qJ(2,3));
t284 = t211 * t221;
t310 = qJ(3,3) * t205;
t148 = t284 + t310;
t139 = 0.1e1 / t148;
t302 = t139 * t211;
t193 = legFrame(3,3);
t165 = sin(t193);
t168 = cos(t193);
t206 = sin(qJ(1,3));
t212 = cos(qJ(1,3));
t112 = -t165 * t206 + t168 * t212;
t175 = t212 * pkin(4);
t133 = t206 * t310 + t175;
t346 = pkin(4) * t206;
t134 = t212 * t310 - t346;
t73 = t112 * t284 - t133 * t165 + t134 * t168;
t113 = t165 * t212 + t168 * t206;
t74 = t113 * t284 + t133 * t168 + t134 * t165;
t58 = (t205 * t218 + (t219 * t74 + t220 * t73) * t302) * t223;
t326 = t221 * t58;
t292 = t205 * t221;
t309 = qJ(3,3) * t211;
t241 = -t292 + t309;
t103 = t148 * t206 + t175;
t106 = t148 * t212 - t346;
t88 = t103 * t168 + t106 * t165;
t322 = t223 * t88;
t85 = -t103 * t165 + t106 * t168;
t323 = t223 * t85;
t61 = -t218 * t223 * t241 + t219 * t322 + t220 * t323;
t356 = t61 - t326;
t207 = sin(qJ(2,2));
t225 = 0.1e1 / qJ(3,2);
t213 = cos(qJ(2,2));
t282 = t213 * t221;
t312 = qJ(3,2) * t207;
t149 = t282 + t312;
t140 = 0.1e1 / t149;
t300 = t140 * t213;
t194 = legFrame(2,3);
t166 = sin(t194);
t169 = cos(t194);
t208 = sin(qJ(1,2));
t214 = cos(qJ(1,2));
t114 = -t166 * t208 + t169 * t214;
t176 = t214 * pkin(4);
t135 = t208 * t312 + t176;
t344 = pkin(4) * t208;
t136 = t214 * t312 - t344;
t75 = t114 * t282 - t135 * t166 + t136 * t169;
t115 = t166 * t214 + t169 * t208;
t76 = t115 * t282 + t135 * t169 + t136 * t166;
t59 = (t207 * t218 + (t219 * t76 + t220 * t75) * t300) * t225;
t325 = t221 * t59;
t289 = t207 * t221;
t311 = qJ(3,2) * t213;
t242 = -t289 + t311;
t104 = t149 * t208 + t176;
t107 = t149 * t214 - t344;
t89 = t104 * t169 + t107 * t166;
t320 = t225 * t89;
t86 = -t104 * t166 + t107 * t169;
t321 = t225 * t86;
t62 = -t218 * t225 * t242 + t219 * t320 + t220 * t321;
t355 = t62 - t325;
t209 = sin(qJ(2,1));
t227 = 0.1e1 / qJ(3,1);
t215 = cos(qJ(2,1));
t280 = t215 * t221;
t314 = qJ(3,1) * t209;
t150 = t280 + t314;
t141 = 0.1e1 / t150;
t298 = t141 * t215;
t195 = legFrame(1,3);
t167 = sin(t195);
t170 = cos(t195);
t210 = sin(qJ(1,1));
t216 = cos(qJ(1,1));
t116 = -t167 * t210 + t170 * t216;
t177 = t216 * pkin(4);
t137 = t210 * t314 + t177;
t342 = pkin(4) * t210;
t138 = t216 * t314 - t342;
t77 = t116 * t280 - t137 * t167 + t138 * t170;
t117 = t167 * t216 + t170 * t210;
t78 = t117 * t280 + t137 * t170 + t138 * t167;
t60 = (t209 * t218 + (t219 * t78 + t220 * t77) * t298) * t227;
t324 = t221 * t60;
t286 = t209 * t221;
t313 = qJ(3,1) * t215;
t243 = -t286 + t313;
t105 = t150 * t210 + t177;
t108 = t150 * t216 - t342;
t90 = t105 * t170 + t108 * t167;
t318 = t227 * t90;
t87 = -t105 * t167 + t108 * t170;
t319 = t227 * t87;
t63 = -t218 * t227 * t243 + t219 * t318 + t220 * t319;
t354 = t63 - t324;
t353 = 0.2e1 * t61;
t352 = 0.2e1 * t62;
t351 = 0.2e1 * t63;
t350 = pkin(4) * t58;
t349 = pkin(4) * t59;
t348 = pkin(4) * t60;
t347 = pkin(4) * t205;
t345 = pkin(4) * t207;
t343 = pkin(4) * t209;
t171 = m(3) * pkin(1) + mrSges(3,1);
t161 = mrSges(2,1) + t171;
t341 = g(3) * t161;
t340 = Ifges(2,1) + Ifges(3,1);
t339 = Ifges(2,6) - Ifges(3,6);
t338 = mrSges(3,2) * t205;
t337 = mrSges(3,2) * t207;
t336 = mrSges(3,2) * t209;
t335 = mrSges(3,3) * qJ(3,1);
t334 = mrSges(3,3) * qJ(3,2);
t333 = mrSges(3,3) * qJ(3,3);
t189 = t211 ^ 2;
t70 = (t112 * t219 - t113 * t220) * t139;
t332 = t189 * t70;
t190 = t213 ^ 2;
t71 = (t114 * t219 - t115 * t220) * t140;
t331 = t190 * t71;
t191 = t215 ^ 2;
t72 = (t116 * t219 - t117 * t220) * t141;
t330 = t191 * t72;
t329 = t211 * t70;
t328 = t213 * t71;
t327 = t215 * t72;
t162 = m(3) * qJ(3,3) + mrSges(3,3);
t317 = t61 * t162;
t163 = m(3) * qJ(3,2) + mrSges(3,3);
t316 = t62 * t163;
t164 = m(3) * qJ(3,1) + mrSges(3,3);
t315 = t63 * t164;
t308 = t112 * t139;
t307 = t113 * t139;
t306 = t114 * t140;
t305 = t115 * t140;
t304 = t116 * t141;
t303 = t117 * t141;
t296 = t171 * t223;
t295 = t171 * t225;
t294 = t171 * t227;
t293 = t205 * t211;
t291 = t205 * t223;
t290 = t207 * t213;
t288 = t207 * t225;
t287 = t209 * t215;
t285 = t209 * t227;
t283 = t211 * t223;
t281 = t213 * t225;
t279 = t215 * t227;
t278 = 2 * mrSges(3,1) * pkin(1);
t277 = -0.2e1 * qJ(3,1) * (t72 * t286 - t348 / 0.2e1);
t276 = -0.2e1 * qJ(3,2) * (t71 * t289 - t349 / 0.2e1);
t275 = -0.2e1 * qJ(3,3) * (t70 * t292 - t350 / 0.2e1);
t274 = Ifges(1,3) + t340;
t273 = pkin(1) * mrSges(3,3) + Ifges(2,4) - Ifges(3,5);
t272 = mrSges(3,2) * t291;
t271 = mrSges(3,2) * t288;
t270 = mrSges(3,2) * t285;
t269 = t73 * t302;
t268 = t74 * t302;
t185 = 0.2e1 * t333;
t222 = qJ(3,3) ^ 2;
t230 = pkin(1) ^ 2;
t247 = Ifges(3,2) + Ifges(2,3) + t278;
t124 = (t222 + t230) * m(3) + t185 + t247;
t97 = (t124 * t205 + t171 * t241) * t223;
t267 = t97 * t302;
t266 = t75 * t300;
t265 = t76 * t300;
t186 = 0.2e1 * t334;
t224 = qJ(3,2) ^ 2;
t125 = (t224 + t230) * m(3) + t186 + t247;
t98 = (t125 * t207 + t171 * t242) * t225;
t264 = t98 * t300;
t263 = t77 * t298;
t262 = t78 * t298;
t187 = 0.2e1 * t335;
t226 = qJ(3,1) ^ 2;
t126 = (t226 + t230) * m(3) + t187 + t247;
t99 = (t126 * t209 + t171 * t243) * t227;
t261 = t99 * t298;
t154 = mrSges(3,2) * qJ(3,3) + t339;
t160 = pkin(1) * mrSges(3,2) - Ifges(3,4) - Ifges(2,5);
t109 = t154 * t211 - t205 * t160;
t260 = t109 * t283;
t155 = mrSges(3,2) * qJ(3,2) + t339;
t110 = t155 * t213 - t207 * t160;
t259 = t110 * t281;
t156 = qJ(3,1) * mrSges(3,2) + t339;
t111 = t156 * t215 - t209 * t160;
t258 = t111 * t279;
t257 = t124 * t283;
t256 = t125 * t281;
t255 = t126 * t279;
t254 = t171 * t283;
t253 = t171 * t281;
t252 = t171 * t279;
t251 = (qJ(3,3) + t221) * (-qJ(3,3) + t221) * t189;
t250 = (qJ(3,2) + t221) * (-qJ(3,2) + t221) * t190;
t249 = (qJ(3,1) + t221) * (-qJ(3,1) + t221) * t191;
t248 = -t230 + (-2 * pkin(1) - pkin(2)) * pkin(2);
t127 = -g(1) * t165 + g(2) * t168;
t130 = g(1) * t168 + g(2) * t165;
t240 = t127 * t206 + t130 * t212;
t128 = -g(1) * t166 + g(2) * t169;
t131 = g(1) * t169 + g(2) * t166;
t239 = t128 * t208 + t131 * t214;
t129 = -g(1) * t167 + g(2) * t170;
t132 = g(1) * t170 + g(2) * t167;
t238 = t129 * t210 + t132 * t216;
t157 = -mrSges(2,2) + t162;
t237 = t157 * t205 + t161 * t211 + mrSges(1,1);
t158 = -mrSges(2,2) + t163;
t236 = t158 * t207 + t161 * t213 + mrSges(1,1);
t159 = -mrSges(2,2) + t164;
t235 = t159 * t209 + t161 * t215 + mrSges(1,1);
t234 = m(3) * t230 + Ifges(2,2) + Ifges(3,3) + t278 - t340;
t228 = pkin(4) ^ 2;
t204 = xDDP(1);
t203 = xDDP(2);
t202 = xDDP(3);
t200 = m(3) * t226;
t198 = m(3) * t224;
t196 = m(3) * t222;
t178 = mrSges(1,2) - mrSges(3,2) - mrSges(2,3);
t174 = t226 + t228;
t173 = t224 + t228;
t172 = t222 + t228;
t144 = qJ(3,1) * t171 + t273;
t143 = qJ(3,2) * t171 + t273;
t142 = qJ(3,3) * t171 + t273;
t120 = -t200 + t234 - 0.2e1 * t335;
t119 = -t198 + t234 - 0.2e1 * t334;
t118 = -t196 + t234 - 0.2e1 * t333;
t102 = (-m(3) * t243 - t171 * t209) * t227;
t101 = (-m(3) * t242 - t171 * t207) * t225;
t100 = (-m(3) * t241 - t171 * t205) * t223;
t96 = (-mrSges(3,2) * t243 + t111) * t285;
t95 = (-mrSges(3,2) * t242 + t110) * t288;
t94 = (-mrSges(3,2) * t241 + t109) * t291;
t93 = t120 * t191 + 0.2e1 * t144 * t287 + t187 + t200 + t274;
t92 = t119 * t190 + 0.2e1 * t143 * t290 + t186 + t198 + t274;
t91 = t118 * t189 + 0.2e1 * t142 * t293 + t185 + t196 + t274;
t69 = t72 ^ 2;
t68 = t71 ^ 2;
t67 = t70 ^ 2;
t66 = t72 * t343;
t65 = t71 * t345;
t64 = t70 * t347;
t57 = t60 ^ 2;
t56 = t59 ^ 2;
t55 = t58 ^ 2;
t48 = m(3) * t319 + (-t117 * t336 - t252 * t77) * t141;
t47 = m(3) * t318 + (t116 * t336 - t252 * t78) * t141;
t46 = m(3) * t321 + (-t115 * t337 - t253 * t75) * t140;
t45 = m(3) * t320 + (t114 * t337 - t253 * t76) * t140;
t44 = m(3) * t323 + (-t113 * t338 - t254 * t73) * t139;
t43 = m(3) * t322 + (t112 * t338 - t254 * t74) * t139;
t42 = t144 * t60;
t41 = t143 * t59;
t40 = t142 * t58;
t39 = -t87 * t294 + (-t111 * t117 + t255 * t77) * t141;
t38 = -t90 * t294 + (t111 * t116 + t255 * t78) * t141;
t37 = -t86 * t295 + (-t110 * t115 + t256 * t75) * t140;
t36 = -t89 * t295 + (t110 * t114 + t256 * t76) * t140;
t35 = -t85 * t296 + (-t109 * t113 + t257 * t73) * t139;
t34 = -t88 * t296 + (t109 * t112 + t257 * t74) * t139;
t33 = t87 * t270 + (-t117 * t93 + t258 * t77) * t141;
t32 = t90 * t270 + (t116 * t93 + t258 * t78) * t141;
t31 = t86 * t271 + (-t115 * t92 + t259 * t75) * t140;
t30 = t89 * t271 + (t114 * t92 + t259 * t76) * t140;
t29 = t85 * t272 + (-t113 * t91 + t260 * t73) * t139;
t28 = t88 * t272 + (t112 * t91 + t260 * t74) * t139;
t27 = t66 + t324;
t26 = t65 + t325;
t25 = t64 + t326;
t21 = t354 * t343;
t20 = t355 * t345;
t19 = t356 * t347;
t18 = (-pkin(4) * t72 + (t243 + t313) * t60 + (t351 - t324) * t209) * t72 * t141;
t17 = (-pkin(4) * t71 + (t242 + t311) * t59 + (t352 - t325) * t207) * t71 * t140;
t16 = (-pkin(4) * t70 + (t241 + t309) * t58 + (t353 - t326) * t205) * t70 * t139;
t15 = (((-t226 + t248) * t60 + t221 * t63) * t60 + t27 * t63 + (t215 * t277 + t21 + (-t174 - t249) * t72 + t243 * t348) * t72) * t227;
t14 = (((-t224 + t248) * t59 + t221 * t62) * t59 + t26 * t62 + (t213 * t276 + t20 + (-t173 - t250) * t71 + t242 * t349) * t71) * t225;
t13 = (((-t222 + t248) * t58 + t221 * t61) * t58 + t25 * t61 + (t211 * t275 + t19 + (-t172 - t251) * t70 + t241 * t350) * t70) * t223;
t12 = ((-t249 * t327 + t191 * t277 + (-t174 * t72 + t21) * t215) * t72 + (-(t66 - t354) * t280 + (pkin(4) * t330 + t209 * t354) * qJ(3,1)) * t60 + (t215 * t27 + t314 * t60) * t63) * t141 * t227;
t11 = ((-t250 * t328 + t190 * t276 + (-t173 * t71 + t20) * t213) * t71 + (-(t65 - t355) * t282 + (pkin(4) * t331 + t207 * t355) * qJ(3,2)) * t59 + (t213 * t26 + t312 * t59) * t62) * t140 * t225;
t10 = ((-t251 * t329 + t189 * t275 + (-t172 * t70 + t19) * t211) * t70 + (-(t64 - t356) * t284 + (pkin(4) * t332 + t205 * t356) * qJ(3,3)) * t58 + (t211 * t25 + t310 * t58) * t61) * t139 * t223;
t9 = t171 * t12 + (t191 * t69 - t57 - t69) * t164 + (t215 * g(3) - t15) * m(3) + (-t69 * t171 * t215 - m(3) * t238 - mrSges(3,2) * t18) * t209;
t8 = t171 * t11 + (t190 * t68 - t56 - t68) * t163 + (g(3) * t213 - t14) * m(3) + (-t68 * t171 * t213 - m(3) * t239 - mrSges(3,2) * t17) * t207;
t7 = t171 * t10 + (t189 * t67 - t55 - t67) * t162 + (t211 * g(3) - t13) * m(3) + (-t67 * t171 * t211 - m(3) * t240 - mrSges(3,2) * t16) * t205;
t6 = -t111 * t18 - t126 * t12 + t171 * t15 + 0.2e1 * t60 * t315 + (-t159 * t238 - t341) * t215 + (-g(3) * t159 + t161 * t238) * t209 + (t120 * t287 - 0.2e1 * t144 * t191 + t144) * t69;
t5 = -t110 * t17 - t125 * t11 + t171 * t14 + 0.2e1 * t59 * t316 + (-t158 * t239 - t341) * t213 + (-g(3) * t158 + t161 * t239) * t207 + (t119 * t290 - 0.2e1 * t143 * t190 + t143) * t68;
t4 = -t109 * t16 - t124 * t10 + t171 * t13 + 0.2e1 * t58 * t317 + (-t157 * t240 - t341) * t211 + (-g(3) * t157 + t161 * t240) * t205 + (t118 * t293 - 0.2e1 * t142 * t189 + t142) * t67;
t3 = -t93 * t18 - t111 * t12 + 0.4e1 * (t42 - t315 / 0.2e1) * t330 + (mrSges(3,2) * t351 - t160 * t60) * t60 * t215 - 0.2e1 * (t42 - t315) * t72 + (-mrSges(3,2) * t15 + 0.2e1 * (-t120 * t60 + t171 * t63) * t327 - t57 * t156) * t209 + (t178 * t216 + t210 * t235) * t132 + (t178 * t210 - t216 * t235) * t129;
t2 = -t92 * t17 - t110 * t11 + 0.4e1 * (t41 - t316 / 0.2e1) * t331 + (mrSges(3,2) * t352 - t160 * t59) * t59 * t213 - 0.2e1 * (t41 - t316) * t71 + (-mrSges(3,2) * t14 + 0.2e1 * (-t119 * t59 + t171 * t62) * t328 - t56 * t155) * t207 + (t178 * t214 + t208 * t236) * t131 + (t178 * t208 - t214 * t236) * t128;
t1 = -t91 * t16 - t109 * t10 + 0.4e1 * (t40 - t317 / 0.2e1) * t332 + (mrSges(3,2) * t353 - t160 * t58) * t58 * t211 - 0.2e1 * (t40 - t317) * t70 + (-mrSges(3,2) * t13 + 0.2e1 * (-t118 * t58 + t171 * t61) * t329 - t55 * t154) * t205 + (t178 * t212 + t206 * t237) * t130 + (t178 * t206 - t237 * t212) * t127;
t22 = [-t1 * t307 - t2 * t305 - t3 * t303 - g(1) * m(4) + (t29 * t308 + t304 * t33 + t306 * t31) * t203 + (-t29 * t307 - t303 * t33 - t305 * t31 + m(4)) * t204 + ((t263 * t39 + t48 * t87) * t204 + (t262 * t39 + t48 * t90) * t203 + (t209 * t39 - t243 * t48) * t202 + t6 * t263 + t87 * t9) * t227 + ((t266 * t37 + t46 * t86) * t204 + (t265 * t37 + t46 * t89) * t203 + (t207 * t37 - t242 * t46) * t202 + t5 * t266 + t86 * t8) * t225 + ((t269 * t35 + t44 * t85) * t204 + (t268 * t35 + t44 * t88) * t203 + (t205 * t35 - t241 * t44) * t202 + t4 * t269 + t85 * t7) * t223; t1 * t308 + t2 * t306 + t3 * t304 - g(2) * m(4) + (-t28 * t307 - t30 * t305 - t303 * t32) * t204 + (t28 * t308 + t30 * t306 + t304 * t32 + m(4)) * t203 + ((t263 * t38 + t47 * t87) * t204 + (t262 * t38 + t47 * t90) * t203 + (t209 * t38 - t243 * t47) * t202 + t6 * t262 + t90 * t9) * t227 + ((t266 * t36 + t45 * t86) * t204 + (t265 * t36 + t45 * t89) * t203 + (t207 * t36 - t242 * t45) * t202 + t5 * t265 + t89 * t8) * t225 + ((t269 * t34 + t43 * t85) * t204 + (t268 * t34 + t43 * t88) * t203 + (t205 * t34 - t241 * t43) * t202 + t4 * t268 + t88 * t7) * t223; (-g(3) + t202) * m(4) + (-t96 * t303 - t95 * t305 - t94 * t307) * t204 + (t96 * t304 + t95 * t306 + t94 * t308) * t203 + ((t102 * t87 + t261 * t77) * t204 + (t102 * t90 + t261 * t78) * t203 + (-t102 * t243 + t209 * t99) * t202 + t209 * t6 - t243 * t9) * t227 + ((t101 * t86 + t264 * t75) * t204 + (t101 * t89 + t264 * t76) * t203 + (-t101 * t242 + t207 * t98) * t202 + t207 * t5 - t242 * t8) * t225 + ((t100 * t85 + t267 * t73) * t204 + (t100 * t88 + t267 * t74) * t203 + (-t100 * t241 + t205 * t97) * t202 + t205 * t4 - t241 * t7) * t223;];
tauX  = t22;
