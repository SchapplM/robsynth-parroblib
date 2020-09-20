% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRRR9V1G3A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:58
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:56:30
% EndTime: 2020-08-06 18:56:36
% DurationCPUTime: 6.52s
% Computational Cost: add. (24843->621), mult. (32331->1033), div. (4665->10), fcn. (26568->36), ass. (0->373)
t205 = cos(pkin(7));
t163 = t205 * pkin(2);
t415 = -0.2e1 * t163;
t413 = 2 * pkin(1);
t414 = t413 / 0.2e1;
t234 = (m(2) + m(3));
t187 = pkin(1) * t234;
t150 = t163 + pkin(1);
t195 = pkin(7) + qJ(3,1);
t162 = cos(t195);
t398 = pkin(3) * t162;
t114 = t150 + t398;
t194 = pkin(7) + qJ(3,2);
t161 = cos(t194);
t399 = pkin(3) * t161;
t113 = t150 + t399;
t193 = pkin(7) + qJ(3,3);
t160 = cos(t193);
t400 = pkin(3) * t160;
t112 = t150 + t400;
t206 = (pkin(5) + qJ(2,3));
t140 = mrSges(3,1) * t206 - Ifges(3,5);
t207 = (pkin(5) + qJ(2,2));
t141 = mrSges(3,1) * t207 - Ifges(3,5);
t208 = (pkin(5) + qJ(2,1));
t142 = mrSges(3,1) * t208 - Ifges(3,5);
t137 = mrSges(3,2) * t206 - Ifges(3,6);
t138 = mrSges(3,2) * t207 - Ifges(3,6);
t139 = mrSges(3,2) * t208 - Ifges(3,6);
t412 = 0.2e1 * pkin(3);
t192 = t205 ^ 2;
t411 = 0.2e1 * t192;
t204 = sin(pkin(7));
t410 = -0.2e1 * t204;
t409 = 0.4e1 * t205;
t210 = (mrSges(2,3) + mrSges(3,3));
t408 = 2 * t210;
t407 = -4 * pkin(5) - 4 * pkin(6);
t406 = pkin(1) * mrSges(3,2);
t154 = 0.1e1 / t160;
t189 = pkin(6) + t206;
t171 = 1 / t189;
t218 = sin(qJ(1,3));
t230 = xDP(3);
t231 = xDP(2);
t232 = xDP(1);
t157 = sin(t193);
t211 = legFrame(3,2);
t177 = cos(t211);
t174 = sin(t211);
t224 = cos(qJ(1,3));
t340 = t174 * t224;
t97 = t157 * t177 - t160 * t340;
t334 = t177 * t224;
t98 = t157 * t174 + t160 * t334;
t67 = (-t218 * t230 + (t231 * t97 + t232 * t98) * t154) * t171;
t65 = pkin(1) * t67;
t158 = sin(t194);
t212 = legFrame(2,2);
t175 = sin(t212);
t178 = cos(t212);
t226 = cos(qJ(1,2));
t332 = t178 * t226;
t100 = t158 * t175 + t161 * t332;
t155 = 0.1e1 / t161;
t190 = pkin(6) + t207;
t172 = 1 / t190;
t220 = sin(qJ(1,2));
t338 = t175 * t226;
t99 = t158 * t178 - t161 * t338;
t68 = (-t220 * t230 + (t100 * t232 + t231 * t99) * t155) * t172;
t66 = pkin(1) * t68;
t405 = pkin(2) * mrSges(3,1);
t404 = pkin(2) * mrSges(3,2);
t159 = sin(t195);
t213 = legFrame(1,2);
t179 = cos(t213);
t176 = sin(t213);
t228 = cos(qJ(1,1));
t336 = t176 * t228;
t101 = t159 * t179 - t162 * t336;
t330 = t179 * t228;
t102 = t159 * t176 + t162 * t330;
t156 = 0.1e1 / t162;
t191 = pkin(6) + t208;
t173 = 1 / t191;
t222 = sin(qJ(1,1));
t69 = (-t222 * t230 + (t101 * t231 + t102 * t232) * t156) * t173;
t64 = t69 * pkin(1);
t227 = cos(qJ(3,1));
t221 = sin(qJ(3,1));
t320 = t221 * t204;
t108 = 0.1e1 / (t227 * t205 - t320);
t348 = t108 * t173;
t247 = -pkin(3) / 0.2e1;
t392 = t227 * pkin(2);
t202 = t227 ^ 2;
t395 = t202 * pkin(3);
t123 = t395 + t392 / 0.2e1 + t247;
t277 = t228 * t320;
t314 = pkin(1) * t228 + t222 * t191;
t261 = pkin(2) * t277 + (t277 * t412 - t314) * t227;
t248 = pkin(2) / 0.2e1;
t341 = (pkin(3) * t227 + t248) * t221;
t164 = pkin(1) * t204;
t344 = (-pkin(3) * t221 + t164) * t227;
t84 = t314 * t320 + (t202 - 0.1e1) * t228 * pkin(3);
t96 = pkin(1) * t221 + (-pkin(3) + t392 + 0.2e1 * t395) * t204;
t50 = (-t123 * t336 + t179 * t341) * t411 + (t261 * t176 + t179 * t96) * t205 + t84 * t176 + t179 * t344;
t44 = t50 * t231 * t348;
t51 = (t123 * t330 + t176 * t341) * t411 + (t176 * t96 - t261 * t179) * t205 - t84 * t179 + t176 * t344;
t45 = t51 * t232 * t348;
t93 = -t114 * t222 + t191 * t228;
t81 = t93 * t173 * t230;
t19 = t64 - 0.2e1 * t45 - 0.2e1 * t44 - 0.2e1 * t81;
t403 = mrSges(3,2) * t19;
t223 = cos(qJ(3,3));
t217 = sin(qJ(3,3));
t322 = t217 * t204;
t106 = 0.1e1 / (t223 * t205 - t322);
t352 = t106 * t171;
t394 = t223 * pkin(2);
t200 = t223 ^ 2;
t397 = t200 * pkin(3);
t121 = t397 + t394 / 0.2e1 + t247;
t279 = t224 * t322;
t316 = pkin(1) * t224 + t218 * t189;
t263 = pkin(2) * t279 + (t279 * t412 - t316) * t223;
t343 = (pkin(3) * t223 + t248) * t217;
t346 = (-pkin(3) * t217 + t164) * t223;
t82 = t316 * t322 + (t200 - 0.1e1) * t224 * pkin(3);
t94 = pkin(1) * t217 + (-pkin(3) + t394 + 0.2e1 * t397) * t204;
t46 = (-t121 * t340 + t177 * t343) * t411 + (t263 * t174 + t177 * t94) * t205 + t82 * t174 + t177 * t346;
t40 = t46 * t231 * t352;
t47 = (t121 * t334 + t174 * t343) * t411 + (t174 * t94 - t263 * t177) * t205 - t82 * t177 + t174 * t346;
t41 = t47 * t232 * t352;
t91 = -t112 * t218 + t189 * t224;
t79 = t91 * t171 * t230;
t20 = t65 - 0.2e1 * t41 - 0.2e1 * t40 - 0.2e1 * t79;
t402 = mrSges(3,2) * t20;
t225 = cos(qJ(3,2));
t219 = sin(qJ(3,2));
t321 = t219 * t204;
t107 = 0.1e1 / (t225 * t205 - t321);
t350 = t107 * t172;
t393 = t225 * pkin(2);
t201 = t225 ^ 2;
t396 = t201 * pkin(3);
t122 = t396 + t393 / 0.2e1 + t247;
t278 = t226 * t321;
t315 = pkin(1) * t226 + t220 * t190;
t262 = pkin(2) * t278 + (t278 * t412 - t315) * t225;
t342 = (pkin(3) * t225 + t248) * t219;
t345 = (-pkin(3) * t219 + t164) * t225;
t83 = t315 * t321 + (t201 - 0.1e1) * t226 * pkin(3);
t95 = pkin(1) * t219 + (-pkin(3) + t393 + 0.2e1 * t396) * t204;
t48 = (-t122 * t338 + t178 * t342) * t411 + (t262 * t175 + t178 * t95) * t205 + t83 * t175 + t178 * t345;
t42 = t48 * t231 * t350;
t49 = (t122 * t332 + t175 * t342) * t411 + (t175 * t95 - t262 * t178) * t205 - t83 * t178 + t175 * t345;
t43 = t49 * t232 * t350;
t92 = -t113 * t220 + t190 * t226;
t80 = t92 * t172 * t230;
t21 = t66 - 0.2e1 * t43 - 0.2e1 * t42 - 0.2e1 * t80;
t401 = mrSges(3,2) * t21;
t209 = Ifges(3,1) - Ifges(3,2);
t391 = 2 * Ifges(2,4) - 2 * Ifges(3,4);
t390 = mrSges(3,1) * t204;
t389 = mrSges(3,2) * t219;
t388 = mrSges(3,2) * t221;
t387 = Ifges(3,4) * t217;
t386 = Ifges(3,4) * t219;
t385 = Ifges(3,4) * t221;
t384 = t106 * t46;
t383 = t106 * t47;
t382 = t107 * t48;
t381 = t107 * t49;
t380 = t108 * t50;
t379 = t108 * t51;
t253 = 0.1e1 / pkin(3);
t88 = (t174 * t232 + t177 * t231) * t253 * t154;
t378 = t88 ^ 2 * t154;
t377 = t154 * t97;
t376 = t154 * t98;
t89 = (t175 * t232 + t178 * t231) * t253 * t155;
t375 = t89 ^ 2 * t155;
t374 = t155 * t99;
t90 = (t176 * t232 + t179 * t231) * t253 * t156;
t373 = t90 ^ 2 * t156;
t372 = t171 * t97;
t371 = t171 * t98;
t370 = t172 * t99;
t369 = t192 * t67;
t368 = t192 * t68;
t367 = t192 * t69;
t366 = t200 * Ifges(3,4);
t365 = t201 * Ifges(3,4);
t364 = t202 * Ifges(3,4);
t363 = t204 * t67;
t362 = t204 * t68;
t361 = t204 * t69;
t165 = t217 * mrSges(3,1);
t360 = t217 * mrSges(3,2);
t166 = t219 * mrSges(3,1);
t167 = t221 * mrSges(3,1);
t359 = m(3) * pkin(2) + mrSges(2,1);
t358 = t100 * t155;
t357 = t100 * t172;
t356 = t101 * t156;
t355 = t101 * t173;
t354 = t102 * t156;
t353 = t102 * t173;
t351 = t106 * t234;
t349 = t107 * t234;
t347 = t108 * t234;
t339 = t174 * t253;
t337 = t175 * t253;
t335 = t176 * t253;
t333 = t177 * t253;
t331 = t178 * t253;
t329 = t179 * t253;
t328 = t209 * t200;
t327 = t209 * t201;
t326 = t209 * t202;
t325 = t209 * t217;
t324 = t209 * t219;
t323 = t209 * t221;
t319 = mrSges(3,2) * t223 + t165;
t318 = mrSges(3,2) * t225 + t166;
t317 = mrSges(3,2) * t227 + t167;
t313 = -2 * t406;
t312 = -0.2e1 * t164;
t310 = mrSges(3,1) * t65;
t309 = mrSges(3,1) * t66;
t308 = mrSges(3,1) * t64;
t307 = 0.4e1 * t387;
t306 = 0.4e1 * t386;
t305 = 0.4e1 * t385;
t303 = mrSges(3,2) * t164;
t302 = pkin(2) * t360;
t301 = pkin(2) * t389;
t300 = pkin(2) * t388;
t299 = -mrSges(1,2) + t210;
t28 = t41 + t40 + t79;
t29 = t43 + t42 + t80;
t30 = t45 + t44 + t81;
t298 = t67 * t366;
t297 = t68 * t365;
t296 = t69 * t364;
t295 = -t406 / 0.2e1;
t294 = -t405 / 0.2e1;
t293 = -t405 / 0.4e1;
t292 = m(3) * pkin(5) + t210;
t272 = mrSges(3,1) * t223 - t360;
t266 = (-t272 - t359) * t205 + (mrSges(2,2) + t319) * t204;
t76 = -t187 + t266;
t291 = t76 * t352;
t271 = mrSges(3,1) * t225 - t389;
t265 = (-t271 - t359) * t205 + (mrSges(2,2) + t318) * t204;
t77 = -t187 + t265;
t290 = t77 * t350;
t270 = mrSges(3,1) * t227 - t388;
t264 = (-t270 - t359) * t205 + (mrSges(2,2) + t317) * t204;
t78 = -t187 + t264;
t289 = t78 * t348;
t288 = t157 * t378;
t287 = t158 * t375;
t286 = t159 * t373;
t73 = (-t137 * t223 - t140 * t217) * t205 + t204 * (t137 * t217 - t140 * t223);
t285 = t218 * t253 * t73;
t74 = (-t138 * t225 - t141 * t219) * t205 + t204 * (t138 * t219 - t141 * t225);
t284 = t220 * t253 * t74;
t75 = (-t139 * t227 - t142 * t221) * t205 + t204 * (t139 * t221 - t142 * t227);
t283 = t222 * t253 * t75;
t282 = (t302 - t209) * t363;
t281 = (t301 - t209) * t362;
t280 = (t300 - t209) * t361;
t276 = t328 * t363;
t275 = t327 * t362;
t274 = t326 * t361;
t254 = pkin(2) ^ 2;
t273 = t254 * m(3) - Ifges(2,1) + Ifges(2,2) + t209;
t118 = g(1) * t177 - g(2) * t174;
t269 = -g(3) * t218 + t118 * t224;
t119 = g(1) * t178 - g(2) * t175;
t268 = -g(3) * t220 + t119 * t226;
t120 = g(1) * t179 - g(2) * t176;
t267 = -g(3) * t222 + t120 * t228;
t260 = m(2) * qJ(2,1) + m(3) * t208 + t299;
t259 = m(2) * qJ(2,2) + m(3) * t207 + t299;
t258 = m(2) * qJ(2,3) + m(3) * t206 + t299;
t255 = pkin(1) ^ 2;
t257 = (2 * mrSges(3,3) * pkin(5)) + (t234 * t255) + Ifges(2,1) + Ifges(3,2) + Ifges(1,3);
t244 = 0.2e1 * pkin(7);
t252 = pkin(3) ^ 2;
t256 = -t254 * cos(t244) - (2 * pkin(6) ^ 2) - t252 - t254 - (2 * t255) + ((-4 * pkin(6) - 2 * pkin(5)) * pkin(5));
t251 = qJ(2,1) ^ 2;
t250 = qJ(2,2) ^ 2;
t249 = qJ(2,3) ^ 2;
t246 = -pkin(5) / 0.2e1;
t245 = -pkin(5) / 0.4e1;
t237 = -Ifges(3,4) / 0.2e1;
t236 = Ifges(3,5) / 0.4e1;
t235 = Ifges(3,6) / 0.2e1;
t216 = xDDP(1);
t215 = xDDP(2);
t214 = xDDP(3);
t198 = mrSges(3,1) * t413;
t197 = 0.2e1 * t405;
t196 = -0.2e1 * t404;
t186 = t405 / 0.4e1;
t185 = -t404 / 0.2e1;
t184 = -t404 / 0.4e1;
t180 = Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1;
t153 = 0.2e1 * t195;
t152 = 0.2e1 * t194;
t151 = 0.2e1 * t193;
t146 = mrSges(1,1) + t187;
t133 = g(3) * t146;
t129 = qJ(2,1) * t234 + t292;
t128 = qJ(2,2) * t234 + t292;
t127 = qJ(2,3) * t234 + t292;
t117 = g(1) * t176 + g(2) * t179;
t116 = g(1) * t175 + g(2) * t178;
t115 = g(1) * t174 + g(2) * t177;
t72 = (-t222 * t78 + t234 * t93) * t173;
t71 = (-t220 * t77 + t234 * t92) * t172;
t70 = (-t218 * t76 + t234 * t91) * t171;
t63 = t69 * t323;
t62 = t68 * t324;
t61 = t67 * t325;
t60 = (Ifges(3,3) * t335 + t75 * t353) * t156;
t59 = (Ifges(3,3) * t329 + t75 * t355) * t156;
t58 = (Ifges(3,3) * t337 + t74 * t357) * t155;
t57 = (Ifges(3,3) * t331 + t74 * t370) * t155;
t56 = (Ifges(3,3) * t339 + t73 * t371) * t154;
t55 = (Ifges(3,3) * t333 + t73 * t372) * t154;
t54 = (-0.2e1 * t326 + (t197 + t305) * t227 - 0.2e1 * t300 + t273) * t192 + (t198 * t227 + (t359 - t388) * t413 + (0.4e1 * t364 + t196 * t227 + 0.2e1 * (t209 * t227 - t405) * t221 + t391) * t204) * t205 + t326 + 0.2e1 * (-t303 - t385) * t227 + (t167 + mrSges(2,2)) * t312 + (t208 ^ 2) * m(3) + m(2) * t251 + qJ(2,1) * t408 + t257;
t53 = (-0.2e1 * t327 + (t197 + t306) * t225 - 0.2e1 * t301 + t273) * t192 + (t198 * t225 + (t359 - t389) * t413 + (0.4e1 * t365 + t196 * t225 + 0.2e1 * (t209 * t225 - t405) * t219 + t391) * t204) * t205 + t327 + 0.2e1 * (-t303 - t386) * t225 + (t166 + mrSges(2,2)) * t312 + (t207 ^ 2) * m(3) + m(2) * t250 + qJ(2,2) * t408 + t257;
t52 = (-0.2e1 * t328 + (t197 + t307) * t223 - 0.2e1 * t302 + t273) * t192 + (t198 * t223 + (t359 - t360) * t413 + (0.4e1 * t366 + t196 * t223 + 0.2e1 * (t209 * t223 - t405) * t217 + t391) * t204) * t205 + t328 + 0.2e1 * (-t303 - t387) * t223 + (t165 + mrSges(2,2)) * t312 + (t206 ^ 2) * m(3) + m(2) * t249 + qJ(2,3) * t408 + t257;
t39 = (-t222 * t54 + t78 * t93) * t173;
t38 = (-t220 * t53 + t77 * t92) * t172;
t37 = (-t218 * t52 + t76 * t91) * t171;
t36 = (t51 * t347 + t78 * t354) * t173;
t35 = (t50 * t347 + t78 * t356) * t173;
t34 = (t49 * t349 + t77 * t358) * t172;
t33 = (t48 * t349 + t77 * t374) * t172;
t32 = (t47 * t351 + t76 * t376) * t171;
t31 = (t46 * t351 + t76 * t377) * t171;
t27 = t51 * t289 + (t75 * t335 + t54 * t353) * t156;
t26 = t50 * t289 + (t75 * t329 + t54 * t355) * t156;
t25 = t49 * t290 + (t74 * t337 + t53 * t357) * t155;
t24 = t48 * t290 + (t74 * t331 + t53 * t370) * t155;
t23 = t47 * t291 + (t73 * t339 + t52 * t371) * t154;
t22 = t46 * t291 + (t73 * t333 + t52 * t372) * t154;
t18 = t66 - t43 / 0.2e1 - t42 / 0.2e1 - t80 / 0.2e1;
t17 = t65 - t41 / 0.2e1 - t40 / 0.2e1 - t79 / 0.2e1;
t16 = t64 - t45 / 0.2e1 - t44 / 0.2e1 - t81 / 0.2e1;
t15 = (-pkin(3) * t373 + (-t64 + 0.2e1 * t30 + (-t163 - t398) * t69) * t69) * t173;
t14 = (-pkin(3) * t375 + (-t66 + 0.2e1 * t29 + (-t163 - t399) * t68) * t68) * t172;
t13 = (-pkin(3) * t378 + (-t65 + 0.2e1 * t28 + (-t163 - t400) * t67) * t67) * t171;
t12 = ((t16 * t415 + (qJ(2,1) * t407 - t252 * cos(t153) - 0.2e1 * t251 + t256) * t69 / 0.2e1 + (t114 + t414) * t30) * t69 + ((t90 * t191 * t159 - 0.2e1 * t16 * t162 + (-cos(t244 + qJ(3,1)) * pkin(2) - t392) * t69) * t69 - (-t69 * t191 * sin(t153) / 0.2e1 + t90 * t114) * t156 * t90) * pkin(3)) * t173;
t11 = ((t18 * t415 + (qJ(2,2) * t407 - t252 * cos(t152) - 0.2e1 * t250 + t256) * t68 / 0.2e1 + (t113 + t414) * t29) * t68 + ((t89 * t190 * t158 - 0.2e1 * t18 * t161 + (-cos(t244 + qJ(3,2)) * pkin(2) - t393) * t68) * t68 - (-t68 * t190 * sin(t152) / 0.2e1 + t89 * t113) * t155 * t89) * pkin(3)) * t172;
t10 = ((t17 * t415 + (qJ(2,3) * t407 - t252 * cos(t151) - 0.2e1 * t249 + t256) * t67 / 0.2e1 + (t112 + t414) * t28) * t67 + ((t88 * t189 * t157 - 0.2e1 * t17 * t160 + (-cos(qJ(3,3) + t244) * pkin(2) - t394) * t67) * t67 - (-t67 * t189 * sin(t151) / 0.2e1 + t88 * t112) * t154 * t88) * pkin(3)) * t171;
t9 = -t75 * t15 + Ifges(3,3) * t286 + (-0.4e1 * (t364 + (t180 * t221 + t184) * t227 + t221 * t293 + t237) * t367 + (-0.2e1 * t274 + (t403 + (t305 + t405) * t361) * t227 - t280 + t19 * t167) * t205 + 0.2e1 * t296 + (t19 * t390 + t63) * t227 - t320 * t403 - Ifges(3,4) * t69) * t69 + (-mrSges(3,1) * t117 + t267 * mrSges(3,2)) * t162 + t159 * (t267 * mrSges(3,1) + mrSges(3,2) * t117);
t8 = -t74 * t14 + Ifges(3,3) * t287 + (-0.4e1 * (t365 + (t180 * t219 + t184) * t225 + t219 * t293 + t237) * t368 + (-0.2e1 * t275 + (t401 + (t306 + t405) * t362) * t225 - t281 + t21 * t166) * t205 + 0.2e1 * t297 + (t21 * t390 + t62) * t225 - t321 * t401 - Ifges(3,4) * t68) * t68 + (-mrSges(3,1) * t116 + t268 * mrSges(3,2)) * t161 + t158 * (t268 * mrSges(3,1) + mrSges(3,2) * t116);
t7 = -t73 * t13 + Ifges(3,3) * t288 + (-0.4e1 * (t366 + (t180 * t217 + t184) * t223 + t217 * t293 + t237) * t369 + (-0.2e1 * t276 + (t402 + (t307 + t405) * t363) * t223 - t282 + t20 * t165) * t205 + 0.2e1 * t298 + (t20 * t390 + t61) * t223 - t322 * t402 - Ifges(3,4) * t67) * t67 + (-mrSges(3,1) * t115 + t269 * mrSges(3,2)) * t160 + t157 * (t269 * mrSges(3,1) + mrSges(3,2) * t115);
t6 = -t78 * t15 - (t129 * t69 + 0.2e1 * (-t270 * t204 - t317 * t205) * t90) * t69 + (-g(3) * t228 - t120 * t222 - t12) * t234;
t5 = -t77 * t14 - (t128 * t68 + 0.2e1 * (-t271 * t204 - t318 * t205) * t89) * t68 + (-g(3) * t226 - t119 * t220 - t11) * t234;
t4 = -t76 * t13 - (t127 * t67 + 0.2e1 * (-t272 * t204 - t319 * t205) * t88) * t67 + (-g(3) * t224 - t118 * t218 - t10) * t234;
t3 = -t54 * t15 - t78 * t12 + t75 * t286 + 0.4e1 * t90 * (0.2e1 * t364 + (t185 + t323) * t227 + t221 * t294 - Ifges(3,4)) * t367 + (t274 + (((-qJ(2,1) / 0.4e1 + t245) * mrSges(3,1) + t236) * t90 + ((t186 + t385) * t410 + t295) * t69) * t227 + t280 / 0.2e1 - t221 * (-t139 * t90 + 0.2e1 * t308) / 0.4e1) * t90 * t409 - 0.4e1 * t90 * t296 - 0.2e1 * ((((-qJ(2,1) / 0.2e1 + t246) * mrSges(3,2) + t235) * t90 + t308) * t204 + t63) * t90 * t227 - (-t142 * t90 + t69 * t313) * t90 * t320 + 0.2e1 * t69 * (Ifges(3,4) * t90 + t129 * t30) + (-t264 * g(3) - t260 * t120 + t133) * t228 - (-t260 * g(3) + (-t146 + t264) * t120) * t222;
t2 = -t53 * t14 - t77 * t11 + t74 * t287 + 0.4e1 * t89 * (0.2e1 * t365 + (t185 + t324) * t225 + t219 * t294 - Ifges(3,4)) * t368 + (t275 + (((-qJ(2,2) / 0.4e1 + t245) * mrSges(3,1) + t236) * t89 + ((t186 + t386) * t410 + t295) * t68) * t225 + t281 / 0.2e1 - t219 * (-t138 * t89 + 0.2e1 * t309) / 0.4e1) * t89 * t409 - 0.4e1 * t89 * t297 - 0.2e1 * ((((-qJ(2,2) / 0.2e1 + t246) * mrSges(3,2) + t235) * t89 + t309) * t204 + t62) * t89 * t225 - (-t141 * t89 + t68 * t313) * t89 * t321 + 0.2e1 * t68 * (Ifges(3,4) * t89 + t128 * t29) + (-t265 * g(3) - t259 * t119 + t133) * t226 - (-t259 * g(3) + (-t146 + t265) * t119) * t220;
t1 = -t52 * t13 - t76 * t10 + t73 * t288 + 0.4e1 * t88 * (0.2e1 * t366 + (t185 + t325) * t223 + t217 * t294 - Ifges(3,4)) * t369 + (t276 + (((-qJ(2,3) / 0.4e1 + t245) * mrSges(3,1) + t236) * t88 + ((t186 + t387) * t410 + t295) * t67) * t223 + t282 / 0.2e1 - t217 * (-t137 * t88 + 0.2e1 * t310) / 0.4e1) * t88 * t409 - 0.4e1 * t88 * t298 - 0.2e1 * ((((-qJ(2,3) / 0.2e1 + t246) * mrSges(3,2) + t235) * t88 + t310) * t204 + t61) * t88 * t223 - (-t140 * t88 + t67 * t313) * t88 * t322 + 0.2e1 * t67 * (Ifges(3,4) * t88 + t127 * t28) + (-t266 * g(3) - t258 * t118 + t133) * t224 - (-t258 * g(3) + (-t146 + t266) * t118) * t218;
t85 = [(-g(1) + t216) * m(4) + ((t27 * t354 + t36 * t379) * t216 + (t27 * t356 + t36 * t380) * t215 + (-t222 * t27 + t36 * t93) * t214 + t3 * t354 + t6 * t379) * t173 + ((t179 * t215 * t60 + (t216 * t60 + t9) * t176) * t156 + (t178 * t215 * t58 + (t216 * t58 + t8) * t175) * t155 + (t177 * t215 * t56 + (t216 * t56 + t7) * t174) * t154) * t253 + ((t25 * t358 + t34 * t381) * t216 + (t25 * t374 + t34 * t382) * t215 + (-t220 * t25 + t34 * t92) * t214 + t2 * t358 + t5 * t381) * t172 + ((t23 * t376 + t32 * t383) * t216 + (t23 * t377 + t32 * t384) * t215 + (-t218 * t23 + t32 * t91) * t214 + t1 * t376 + t4 * t383) * t171; (-g(2) + t215) * m(4) + ((t26 * t354 + t35 * t379) * t216 + (t26 * t356 + t35 * t380) * t215 + (-t222 * t26 + t35 * t93) * t214 + t3 * t356 + t6 * t380) * t173 + ((t176 * t216 * t59 + (t215 * t59 + t9) * t179) * t156 + (t175 * t216 * t57 + (t215 * t57 + t8) * t178) * t155 + (t174 * t216 * t55 + (t215 * t55 + t7) * t177) * t154) * t253 + ((t24 * t358 + t33 * t381) * t216 + (t24 * t374 + t33 * t382) * t215 + (-t220 * t24 + t33 * t92) * t214 + t2 * t374 + t5 * t382) * t172 + ((t22 * t376 + t31 * t383) * t216 + (t22 * t377 + t31 * t384) * t215 + (-t218 * t22 + t31 * t91) * t214 + t1 * t377 + t4 * t384) * t171; (-g(3) + t214) * m(4) + ((-t222 * t39 + t72 * t93) * t214 - t222 * t3 + t93 * t6 + (t215 * t50 + t216 * t51) * t72 * t108 + ((t102 * t39 - t176 * t283) * t216 + (t101 * t39 - t179 * t283) * t215) * t156) * t173 + ((-t220 * t38 + t71 * t92) * t214 - t220 * t2 + t92 * t5 + (t215 * t48 + t216 * t49) * t71 * t107 + ((t100 * t38 - t175 * t284) * t216 + (-t178 * t284 + t38 * t99) * t215) * t155) * t172 + ((-t218 * t37 + t70 * t91) * t214 - t218 * t1 + t91 * t4 + (t215 * t46 + t216 * t47) * t70 * t106 + ((-t174 * t285 + t37 * t98) * t216 + (-t177 * t285 + t37 * t97) * t215) * t154) * t171;];
tauX  = t85;
