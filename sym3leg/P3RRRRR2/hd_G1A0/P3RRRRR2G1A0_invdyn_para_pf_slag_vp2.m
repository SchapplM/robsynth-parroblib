% Calculate vector of inverse dynamics forces for parallel robot
% P3RRRRR2G1P1A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2020-03-09 21:05
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:04:12
% EndTime: 2020-03-09 21:04:17
% DurationCPUTime: 4.75s
% Computational Cost: add. (15498->454), mult. (12207->770), div. (7578->18), fcn. (11280->72), ass. (0->360)
t394 = 2 * pkin(1);
t393 = -2 * pkin(1);
t392 = g(3) * mrSges(3,1);
t391 = g(3) * mrSges(3,2);
t201 = sin(qJ(2,3));
t390 = pkin(1) * t201;
t203 = sin(qJ(2,2));
t389 = pkin(1) * t203;
t205 = sin(qJ(2,1));
t388 = pkin(1) * t205;
t207 = cos(qJ(2,3));
t387 = pkin(1) * t207;
t209 = cos(qJ(2,2));
t386 = pkin(1) * t209;
t211 = cos(qJ(2,1));
t385 = pkin(1) * t211;
t206 = cos(qJ(3,3));
t176 = t206 ^ 2;
t384 = pkin(2) * t176;
t208 = cos(qJ(3,2));
t180 = t208 ^ 2;
t383 = pkin(2) * t180;
t210 = cos(qJ(3,1));
t184 = t210 ^ 2;
t382 = pkin(2) * t184;
t381 = Ifges(3,1) + Ifges(2,3);
t192 = legFrame(3,3);
t161 = sin(t192);
t164 = cos(t192);
t100 = -g(1) * t161 + g(2) * t164;
t380 = mrSges(3,1) * t100;
t193 = legFrame(2,3);
t162 = sin(t193);
t165 = cos(t193);
t101 = -g(1) * t162 + g(2) * t165;
t379 = mrSges(3,1) * t101;
t194 = legFrame(1,3);
t163 = sin(t194);
t166 = cos(t194);
t102 = -g(1) * t163 + g(2) * t166;
t378 = mrSges(3,1) * t102;
t377 = mrSges(3,2) * t100;
t376 = mrSges(3,2) * t101;
t375 = mrSges(3,2) * t102;
t103 = g(1) * t164 + g(2) * t161;
t374 = mrSges(3,2) * t103;
t104 = g(1) * t165 + g(2) * t162;
t373 = mrSges(3,2) * t104;
t105 = g(1) * t166 + g(2) * t163;
t372 = mrSges(3,2) * t105;
t200 = sin(qJ(3,3));
t371 = mrSges(3,2) * t200;
t202 = sin(qJ(3,2));
t370 = mrSges(3,2) * t202;
t204 = sin(qJ(3,1));
t369 = mrSges(3,2) * t204;
t368 = Ifges(3,4) * t200;
t367 = Ifges(3,4) * t202;
t366 = Ifges(3,4) * t204;
t178 = 0.1e1 / t206 ^ 2;
t365 = Ifges(3,6) * t178;
t182 = 0.1e1 / t208 ^ 2;
t364 = Ifges(3,6) * t182;
t186 = 0.1e1 / t210 ^ 2;
t363 = Ifges(3,6) * t186;
t212 = xDP(3);
t188 = t212 ^ 2;
t333 = t188 / pkin(2) ^ 2;
t265 = t333 / 0.2e1;
t220 = 1 / pkin(1);
t217 = 0.1e1 / pkin(2);
t339 = t178 * t217;
t286 = (pkin(2) * t206 + t387) * t339;
t173 = 0.1e1 / t201;
t346 = t173 * t200;
t239 = t286 * t346;
t229 = t220 * t239;
t226 = t212 * t229;
t213 = xDP(2);
t310 = qJ(2,3) + qJ(3,3);
t311 = qJ(2,3) - qJ(3,3);
t109 = 0.1e1 / (sin(t310) + sin(t311));
t316 = t217 * t220;
t291 = t109 * t316;
t152 = qJ(1,3) + t192;
t145 = qJ(2,3) + t152;
t136 = qJ(3,3) + t145;
t137 = -qJ(3,3) + t145;
t76 = sin(t152) * t393 + (-sin(t137) - sin(t136)) * pkin(2);
t263 = t76 * t291;
t61 = t213 * t263;
t214 = xDP(1);
t79 = cos(t152) * t393 + (-cos(t137) - cos(t136)) * pkin(2);
t262 = t79 * t291;
t64 = t214 * t262;
t46 = t61 + t64 - t226;
t177 = 0.1e1 / t206;
t345 = t173 * t220;
t277 = t200 * t345;
t257 = t177 * t277;
t133 = cos(t145);
t280 = t133 * t345;
t130 = sin(t145);
t283 = t130 * t345;
t70 = t212 * t257 + t213 * t283 + t214 * t280;
t361 = (t178 * t265 + (t70 + t46 / 0.2e1) * t46) * t201;
t337 = t182 * t217;
t285 = (pkin(2) * t208 + t386) * t337;
t174 = 0.1e1 / t203;
t344 = t174 * t202;
t238 = t285 * t344;
t228 = t220 * t238;
t225 = t212 * t228;
t312 = qJ(2,2) + qJ(3,2);
t313 = qJ(2,2) - qJ(3,2);
t110 = 0.1e1 / (sin(t312) + sin(t313));
t289 = t110 * t316;
t169 = qJ(1,2) + t312;
t138 = t193 + t169;
t170 = qJ(1,2) + t313;
t139 = t193 + t170;
t153 = qJ(1,2) + t193;
t77 = sin(t153) * t393 + (-sin(t139) - sin(t138)) * pkin(2);
t261 = t77 * t289;
t62 = t213 * t261;
t80 = cos(t153) * t393 + (-cos(t139) - cos(t138)) * pkin(2);
t260 = t80 * t289;
t65 = t214 * t260;
t47 = t62 + t65 - t225;
t181 = 0.1e1 / t208;
t343 = t174 * t220;
t276 = t202 * t343;
t256 = t181 * t276;
t190 = qJ(2,2) + qJ(1,2);
t146 = t193 + t190;
t134 = cos(t146);
t279 = t134 * t343;
t131 = sin(t146);
t282 = t131 * t343;
t71 = t212 * t256 + t213 * t282 + t214 * t279;
t360 = (t182 * t265 + (t71 + t47 / 0.2e1) * t47) * t203;
t335 = t186 * t217;
t284 = (pkin(2) * t210 + t385) * t335;
t175 = 0.1e1 / t205;
t342 = t175 * t204;
t237 = t284 * t342;
t227 = t220 * t237;
t224 = t212 * t227;
t314 = qJ(2,1) + qJ(3,1);
t315 = qJ(2,1) - qJ(3,1);
t111 = 0.1e1 / (sin(t314) + sin(t315));
t287 = t111 * t316;
t154 = qJ(1,1) + t194;
t147 = qJ(2,1) + t154;
t140 = qJ(3,1) + t147;
t141 = -qJ(3,1) + t147;
t78 = sin(t154) * t393 + (-sin(t141) - sin(t140)) * pkin(2);
t259 = t78 * t287;
t63 = t213 * t259;
t81 = cos(t154) * t393 + (-cos(t141) - cos(t140)) * pkin(2);
t258 = t81 * t287;
t66 = t214 * t258;
t48 = t63 + t66 - t224;
t185 = 0.1e1 / t210;
t341 = t175 * t220;
t275 = t204 * t341;
t255 = t185 * t275;
t135 = cos(t147);
t278 = t135 * t341;
t132 = sin(t147);
t281 = t132 * t341;
t72 = t212 * t255 + t213 * t281 + t214 * t278;
t359 = (t186 * t265 + (t72 + t48 / 0.2e1) * t48) * t205;
t34 = t46 + t70;
t358 = t206 * t34;
t35 = t47 + t71;
t357 = t208 * t35;
t36 = t48 + t72;
t356 = t210 * t36;
t355 = t109 * t217;
t354 = t110 * t217;
t353 = t111 * t217;
t352 = t130 * t173;
t351 = t131 * t174;
t350 = t132 * t175;
t349 = t133 * t173;
t348 = t134 * t174;
t347 = t135 * t175;
t340 = t177 * t217;
t338 = t181 * t217;
t336 = t185 * t217;
t334 = t188 * t217;
t195 = Ifges(3,2) - Ifges(3,1);
t332 = t195 * t176;
t331 = t195 * t180;
t330 = t195 * t184;
t189 = qJ(1,3) + qJ(2,3);
t158 = cos(t189);
t196 = mrSges(3,3) - mrSges(2,2);
t329 = t196 * t158;
t159 = cos(t190);
t328 = t196 * t159;
t191 = qJ(1,1) + qJ(2,1);
t160 = cos(t191);
t327 = t196 * t160;
t326 = t196 * t207;
t325 = t196 * t209;
t324 = t196 * t211;
t197 = xDDP(3);
t323 = t197 * t220;
t198 = xDDP(2);
t322 = t198 * t220;
t199 = xDDP(1);
t321 = t199 * t220;
t320 = t200 * t201;
t319 = t202 * t203;
t318 = t204 * t205;
t317 = t212 * t217;
t309 = 0.2e1 * t368;
t308 = 0.2e1 * t367;
t307 = 0.2e1 * t366;
t149 = mrSges(3,1) * t387;
t150 = mrSges(3,1) * t386;
t151 = mrSges(3,1) * t385;
t306 = t207 * t384;
t305 = t209 * t383;
t304 = t211 * t382;
t303 = Ifges(3,5) * t333;
t302 = Ifges(3,3) * t333;
t301 = t76 * t355;
t300 = t79 * t355;
t251 = t332 + t381;
t85 = t206 * t309 + t251;
t299 = t85 * t355;
t298 = t77 * t354;
t297 = t80 * t354;
t250 = t331 + t381;
t86 = t208 * t308 + t250;
t296 = t86 * t354;
t295 = t78 * t353;
t294 = t81 * t353;
t249 = t330 + t381;
t87 = t210 * t307 + t249;
t293 = t87 * t353;
t118 = Ifges(3,5) * t200 + Ifges(3,6) * t206;
t292 = t118 * t355;
t119 = Ifges(3,5) * t202 + Ifges(3,6) * t208;
t290 = t119 * t354;
t120 = Ifges(3,5) * t204 + Ifges(3,6) * t210;
t288 = t120 * t353;
t274 = t177 * t317;
t273 = t181 * t317;
t272 = t185 * t317;
t271 = t200 * t333;
t270 = t202 * t333;
t269 = t204 * t333;
t268 = t212 * t320;
t267 = t212 * t319;
t266 = t212 * t318;
t264 = 0.4e1 * Ifges(3,4) * t317;
t254 = t207 * t274;
t253 = t209 * t273;
t252 = t211 * t272;
t248 = -0.2e1 * t34 * t274;
t247 = -0.2e1 * t35 * t273;
t246 = -0.2e1 * t36 * t272;
t245 = t34 * t254;
t244 = t35 * t253;
t243 = t36 * t252;
t242 = mrSges(3,1) * t206 - t371;
t241 = mrSges(3,1) * t208 - t370;
t240 = mrSges(3,1) * t210 - t369;
t215 = (m(2) + m(3));
t219 = pkin(1) ^ 2;
t236 = (t215 * t219) + Ifges(1,3) + t381;
t235 = t264 * t358 + Ifges(3,4) * t248 + (t195 * t200 * t248 + t178 * t303) * t206;
t234 = t264 * t357 + Ifges(3,4) * t247 + (t195 * t202 * t247 + t182 * t303) * t208;
t233 = t264 * t356 + Ifges(3,4) * t246 + (t195 * t204 * t246 + t186 * t303) * t210;
t232 = mrSges(2,1) + t242;
t231 = mrSges(2,1) + t241;
t230 = mrSges(2,1) + t240;
t223 = t210 * t184;
t222 = t208 * t180;
t221 = t206 * t176;
t216 = pkin(2) ^ 2;
t187 = 0.1e1 / t223;
t183 = 0.1e1 / t222;
t179 = 0.1e1 / t221;
t172 = qJ(1,1) + t315;
t171 = qJ(1,1) + t314;
t168 = qJ(1,3) + t311;
t167 = qJ(1,3) + t310;
t157 = sin(t191);
t156 = sin(t190);
t155 = sin(t189);
t148 = pkin(1) * t215 + mrSges(1,1);
t144 = t196 * t388;
t143 = t196 * t389;
t142 = t196 * t390;
t114 = (-mrSges(2,1) + t369) * t385;
t113 = (-mrSges(2,1) + t370) * t386;
t112 = (-mrSges(2,1) + t371) * t387;
t99 = mrSges(2,1) * t105;
t98 = mrSges(2,1) * t104;
t97 = mrSges(2,1) * t103;
t96 = mrSges(3,1) * t105;
t95 = mrSges(3,1) * t104;
t94 = mrSges(3,1) * t103;
t84 = t210 * (-mrSges(3,2) * t388 + Ifges(3,6)) - t204 * (mrSges(3,1) * t388 - Ifges(3,5));
t83 = t208 * (-mrSges(3,2) * t389 + Ifges(3,6)) - t202 * (mrSges(3,1) * t389 - Ifges(3,5));
t82 = t206 * (-mrSges(3,2) * t390 + Ifges(3,6)) - t200 * (mrSges(3,1) * t390 - Ifges(3,5));
t75 = (t151 + t307) * t210 - t114 + t144 + t249;
t74 = (t150 + t308) * t208 - t113 + t143 + t250;
t73 = (t149 + t309) * t206 - t112 + t142 + t251;
t69 = t72 ^ 2;
t68 = t71 ^ 2;
t67 = t70 ^ 2;
t60 = t330 + 0.2e1 * (t151 + t366) * t210 - 0.2e1 * t114 + 0.2e1 * t144 + t236;
t59 = t331 + 0.2e1 * (t150 + t367) * t208 - 0.2e1 * t113 + 0.2e1 * t143 + t236;
t58 = t332 + 0.2e1 * (t149 + t368) * t206 - 0.2e1 * t112 + 0.2e1 * t142 + t236;
t57 = t120 * t336 + (t185 * t75 - t87 * t284) * t275;
t56 = t119 * t338 + (t181 * t74 - t86 * t285) * t276;
t55 = t118 * t340 + (t177 * t73 - t85 * t286) * t277;
t54 = (t81 * t293 + t75 * t347) * t220;
t53 = (t80 * t296 + t74 * t348) * t220;
t52 = (t79 * t299 + t73 * t349) * t220;
t51 = (t78 * t293 + t75 * t350) * t220;
t50 = (t77 * t296 + t74 * t351) * t220;
t49 = (t76 * t299 + t73 * t352) * t220;
t45 = (t75 * t294 + t60 * t347) * t220;
t44 = (t74 * t297 + t59 * t348) * t220;
t43 = (t73 * t300 + t58 * t349) * t220;
t42 = (t75 * t295 + t60 * t350) * t220;
t41 = (t74 * t298 + t59 * t351) * t220;
t40 = (t73 * t301 + t58 * t352) * t220;
t39 = t84 * t336 + (t185 * t60 - t75 * t284) * t275;
t38 = t83 * t338 + (t181 * t59 - t74 * t285) * t276;
t37 = t82 * t340 + (t177 * t58 - t73 * t286) * t277;
t33 = t36 ^ 2;
t32 = t35 ^ 2;
t31 = t34 ^ 2;
t30 = t66 / 0.2e1 + t63 / 0.2e1 - t224 / 0.2e1 + t72;
t29 = t65 / 0.2e1 + t62 / 0.2e1 - t225 / 0.2e1 + t71;
t28 = t64 / 0.2e1 + t61 / 0.2e1 - t226 / 0.2e1 + t70;
t27 = t216 * t36 * t223;
t26 = t216 * t35 * t222;
t25 = t216 * t34 * t221;
t12 = ((-t210 * t72 * t385 - t36 * t382) * t185 * t72 - pkin(2) * t48 * t356 - t187 * t334) * t341;
t11 = ((-t208 * t71 * t386 - t35 * t383) * t181 * t71 - pkin(2) * t47 * t357 - t183 * t334) * t343;
t10 = ((-t206 * t70 * t387 - t34 * t384) * t177 * t70 - pkin(2) * t46 * t358 - t179 * t334) * t345;
t9 = (((-pkin(1) * t36 * t318 + t185 * t212) * t210 + pkin(1) * t252) * t187 * t317 + ((t27 + t30 * t304 * t394 + (-pkin(1) * t185 * t266 + t219 * t72) * t210) * t72 + (t27 + (t36 * t304 - t266) * pkin(1)) * t48) * t335) * t341;
t8 = (((-pkin(1) * t35 * t319 + t181 * t212) * t208 + pkin(1) * t253) * t183 * t317 + ((t26 + t29 * t305 * t394 + (-pkin(1) * t181 * t267 + t219 * t71) * t208) * t71 + (t26 + (t35 * t305 - t267) * pkin(1)) * t47) * t337) * t343;
t7 = (((-pkin(1) * t34 * t320 + t177 * t212) * t206 + pkin(1) * t254) * t179 * t317 + ((t25 + t28 * t306 * t394 + (-pkin(1) * t177 * t268 + t219 * t70) * t206) * t70 + (t25 + (t34 * t306 - t268) * pkin(1)) * t46) * t339) * t345;
t6 = -t75 * t12 - t87 * t9 + t99 * t157 + (t120 * t187 - t363) * t269 + (t240 * t157 - t327) * t105 + (-t196 * t157 - t230 * t160) * t102 + (t230 * t205 - t324) * t69 * pkin(1) + t233;
t5 = -t74 * t11 - t86 * t8 + t98 * t156 + (t119 * t183 - t364) * t270 + (t241 * t156 - t328) * t104 + (-t196 * t156 - t231 * t159) * t101 + (t231 * t203 - t325) * t68 * pkin(1) + t234;
t4 = -t73 * t10 - t85 * t7 + t97 * t155 + (t118 * t179 - t365) * t271 + (t242 * t155 - t329) * t103 + (-t196 * t155 - t232 * t158) * t100 + (t232 * t201 - t326) * t67 * pkin(1) + t235;
t3 = -t60 * t12 - t75 * t9 + (-t372 - t378) * cos(t172) / 0.2e1 + (t96 - t375) * sin(t172) / 0.2e1 + (t372 - t378) * cos(t171) / 0.2e1 + (t96 + t375) * sin(t171) / 0.2e1 - mrSges(2,1) * t102 * t160 + (mrSges(1,2) * t102 + t105 * t148) * sin(qJ(1,1)) + (t187 * t84 - t363) * t269 - t105 * t327 + (-t102 * t196 + t99) * t157 + (mrSges(1,2) * t105 - t102 * t148) * cos(qJ(1,1)) + ((-mrSges(3,1) * t359 - mrSges(3,2) * t243) * t210 + (-mrSges(3,1) * t243 + mrSges(3,2) * t359) * t204 + (-mrSges(2,1) * t205 + t324) * t48 * t30) * t394 + t233;
t2 = -t59 * t11 - t74 * t8 + (-t373 - t379) * cos(t170) / 0.2e1 + (t95 - t376) * sin(t170) / 0.2e1 + (t373 - t379) * cos(t169) / 0.2e1 + (t95 + t376) * sin(t169) / 0.2e1 - mrSges(2,1) * t101 * t159 + (mrSges(1,2) * t101 + t104 * t148) * sin(qJ(1,2)) + (t183 * t83 - t364) * t270 - t104 * t328 + (-t101 * t196 + t98) * t156 + (mrSges(1,2) * t104 - t101 * t148) * cos(qJ(1,2)) + ((-mrSges(3,1) * t360 - mrSges(3,2) * t244) * t208 + (-mrSges(3,1) * t244 + mrSges(3,2) * t360) * t202 + (-mrSges(2,1) * t203 + t325) * t47 * t29) * t394 + t234;
t1 = -t58 * t10 - t73 * t7 + (-t374 - t380) * cos(t168) / 0.2e1 + (t94 - t377) * sin(t168) / 0.2e1 + (t374 - t380) * cos(t167) / 0.2e1 + (t94 + t377) * sin(t167) / 0.2e1 - mrSges(2,1) * t100 * t158 + (mrSges(1,2) * t100 + t103 * t148) * sin(qJ(1,3)) + (t179 * t82 - t365) * t271 - t103 * t329 + (-t100 * t196 + t97) * t155 + (mrSges(1,2) * t103 - t100 * t148) * cos(qJ(1,3)) + ((-mrSges(3,1) * t361 - mrSges(3,2) * t245) * t206 + (-mrSges(3,1) * t245 + mrSges(3,2) * t361) * t200 + (-mrSges(2,1) * t201 + t326) * t46 * t28) * t394 + t235;
t13 = [t1 * t280 + t2 * t279 + t6 * t258 + t5 * t260 + t4 * t262 + t3 * t278 + (-g(1) + t199) * m(4) + (-t54 * t237 + (t45 * t342 + (t81 * t288 + t84 * t347) * t217) * t185 - t53 * t238 + (t44 * t344 + (t80 * t290 + t83 * t348) * t217) * t181 - t52 * t239 + (t43 * t346 + (t79 * t292 + t82 * t349) * t217) * t177) * t323 + (t54 * t295 + t53 * t298 + t52 * t301 + t45 * t350 + t44 * t351 + t43 * t352) * t322 + (t54 * t294 + t53 * t297 + t52 * t300 + t45 * t347 + t44 * t348 + t43 * t349) * t321; t1 * t283 + t2 * t282 + t6 * t259 + t5 * t261 + t4 * t263 + t3 * t281 + (-g(2) + t198) * m(4) + (-t51 * t237 + (t42 * t342 + (t78 * t288 + t84 * t350) * t217) * t185 - t50 * t238 + (t41 * t344 + (t77 * t290 + t83 * t351) * t217) * t181 - t49 * t239 + (t40 * t346 + (t76 * t292 + t82 * t352) * t217) * t177) * t323 + (t51 * t295 + t50 * t298 + t49 * t301 + t42 * t350 + t41 * t351 + t40 * t352) * t322 + (t51 * t294 + t50 * t297 + t49 * t300 + t42 * t347 + t41 * t348 + t40 * t349) * t321; t3 * t255 - t6 * t227 + (-t84 * t12 - t120 * t9 + (mrSges(3,2) * t69 * t385 - t392) * t210 + (t195 * t33 * t210 + t69 * t151 + t187 * t302 + t391) * t204 + (t102 * t157 + t105 * t160) * (mrSges(3,1) * t204 + mrSges(3,2) * t210) + (-0.2e1 * t184 + 0.1e1) * t33 * Ifges(3,4)) * t336 + t2 * t256 - t5 * t228 + (-t83 * t11 - t119 * t8 + (mrSges(3,2) * t68 * t386 - t392) * t208 + (t195 * t32 * t208 + t68 * t150 + t183 * t302 + t391) * t202 + (t101 * t156 + t104 * t159) * (mrSges(3,1) * t202 + mrSges(3,2) * t208) + (-0.2e1 * t180 + 0.1e1) * t32 * Ifges(3,4)) * t338 + t1 * t257 - t4 * t229 + (-t82 * t10 - t118 * t7 + (mrSges(3,2) * t67 * t387 - t392) * t206 + (t195 * t31 * t206 + t67 * t149 + t179 * t302 + t391) * t200 + (t100 * t155 + t103 * t158) * (mrSges(3,1) * t200 + mrSges(3,2) * t206) + (-0.2e1 * t176 + 0.1e1) * t31 * Ifges(3,4)) * t340 - g(3) * m(4) + (t57 * t295 + t56 * t298 + t55 * t301 + t39 * t350 + t38 * t351 + t37 * t352) * t322 + (t57 * t294 + t56 * t297 + t55 * t300 + t39 * t347 + t38 * t348 + t37 * t349) * t321 + ((t39 * t185 - t57 * t284 + (-t120 * t284 + t185 * t84) * t336) * t275 + (t38 * t181 - t56 * t285 + (-t119 * t285 + t181 * t83) * t338) * t276 + (t37 * t177 - t55 * t286 + (-t118 * t286 + t177 * t82) * t340) * t277 + m(4) + (t177 ^ 2 + t181 ^ 2 + t185 ^ 2) * Ifges(3,3) * t217 ^ 2) * t197;];
tauX  = t13;
