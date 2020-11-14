% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR8V2G3A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2020-08-06 18:05
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:04:25
% EndTime: 2020-08-06 18:04:34
% DurationCPUTime: 8.62s
% Computational Cost: add. (39108->481), mult. (87021->870), div. (5112->7), fcn. (84471->22), ass. (0->318)
t195 = sin(qJ(2,1));
t201 = cos(qJ(2,1));
t207 = pkin(7) + pkin(6);
t135 = pkin(2) * t195 - t201 * t207;
t180 = sin(pkin(4));
t182 = cos(pkin(4));
t194 = sin(qJ(3,1));
t284 = t182 * t194;
t105 = pkin(3) * t284 + t135 * t180;
t200 = cos(qJ(3,1));
t296 = t180 * t195;
t178 = t200 ^ 2;
t358 = pkin(3) * t178;
t84 = 0.1e1 / (pkin(2) * t284 + t105 * t200 + t296 * t358);
t193 = sin(qJ(2,2));
t199 = cos(qJ(2,2));
t134 = pkin(2) * t193 - t199 * t207;
t192 = sin(qJ(3,2));
t286 = t182 * t192;
t104 = pkin(3) * t286 + t134 * t180;
t198 = cos(qJ(3,2));
t298 = t180 * t193;
t177 = t198 ^ 2;
t359 = pkin(3) * t177;
t83 = 0.1e1 / (pkin(2) * t286 + t104 * t198 + t298 * t359);
t191 = sin(qJ(2,3));
t197 = cos(qJ(2,3));
t133 = pkin(2) * t191 - t197 * t207;
t190 = sin(qJ(3,3));
t288 = t182 * t190;
t103 = pkin(3) * t288 + t133 * t180;
t196 = cos(qJ(3,3));
t300 = t180 * t191;
t176 = t196 ^ 2;
t360 = pkin(3) * t176;
t82 = 0.1e1 / (pkin(2) * t288 + t103 * t196 + t300 * t360);
t186 = legFrame(1,2);
t166 = sin(t186);
t169 = cos(t186);
t129 = g(1) * t166 + g(2) * t169;
t132 = g(1) * t169 - g(2) * t166;
t181 = cos(pkin(8));
t158 = g(3) * t181;
t179 = sin(pkin(8));
t239 = -t132 * t179 - t158;
t381 = t129 * t180 + t182 * t239;
t185 = legFrame(2,2);
t165 = sin(t185);
t168 = cos(t185);
t128 = g(1) * t165 + g(2) * t168;
t131 = g(1) * t168 - g(2) * t165;
t241 = -t131 * t179 - t158;
t380 = t128 * t180 + t182 * t241;
t184 = legFrame(3,2);
t164 = sin(t184);
t167 = cos(t184);
t127 = g(1) * t164 + g(2) * t167;
t130 = g(1) * t167 - g(2) * t164;
t243 = -t130 * t179 - t158;
t379 = t127 * t180 + t182 * t243;
t204 = xDP(3);
t211 = 0.1e1 / pkin(3);
t205 = xDP(2);
t206 = xDP(1);
t234 = t164 * t205 - t167 * t206;
t136 = pkin(2) * t197 + t191 * t207;
t274 = t197 * t181;
t289 = t181 * t182;
t357 = pkin(3) * t196;
t76 = (t179 * t191 - t182 * t274) * t357 - t136 * t289 + t133 * t179;
t275 = t197 * t179;
t302 = t179 * t182;
t79 = (t181 * t191 + t182 * t275) * t357 + t136 * t302 + t133 * t181;
t58 = (t204 * t79 - t234 * t76) * t82 * t211;
t378 = -0.2e1 * t58;
t233 = t165 * t205 - t168 * t206;
t137 = pkin(2) * t199 + t193 * t207;
t272 = t199 * t181;
t356 = pkin(3) * t198;
t77 = (t179 * t193 - t182 * t272) * t356 - t137 * t289 + t134 * t179;
t273 = t199 * t179;
t80 = (t181 * t193 + t182 * t273) * t356 + t137 * t302 + t134 * t181;
t59 = (t204 * t80 - t233 * t77) * t83 * t211;
t377 = -0.2e1 * t59;
t232 = t166 * t205 - t169 * t206;
t138 = pkin(2) * t201 + t195 * t207;
t270 = t201 * t181;
t355 = pkin(3) * t200;
t78 = (t179 * t195 - t182 * t270) * t355 - t138 * t289 + t135 * t179;
t271 = t201 * t179;
t81 = (t181 * t195 + t182 * t271) * t355 + t138 * t302 + t135 * t181;
t60 = (t204 * t81 - t232 * t78) * t84 * t211;
t376 = -0.2e1 * t60;
t249 = t196 * mrSges(3,1) - mrSges(3,2) * t190;
t248 = t198 * mrSges(3,1) - mrSges(3,2) * t192;
t247 = t200 * mrSges(3,1) - mrSges(3,2) * t194;
t183 = Ifges(3,1) - Ifges(3,2);
t202 = mrSges(3,2) * pkin(2);
t366 = -2 * Ifges(3,4);
t369 = Ifges(3,4) + t178 * t366 + (-t183 * t194 + t202) * t200;
t368 = Ifges(3,4) + t177 * t366 + (-t183 * t192 + t202) * t198;
t367 = Ifges(3,4) + t176 * t366 + (-t183 * t190 + t202) * t196;
t203 = mrSges(3,1) * pkin(2);
t365 = pkin(3) * t58;
t364 = pkin(3) * t59;
t363 = pkin(3) * t60;
t159 = pkin(6) * mrSges(3,2) - Ifges(3,6);
t362 = -t159 / 0.2e1;
t160 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t361 = t160 / 0.2e1;
t354 = t179 * g(3);
t212 = pkin(2) ^ 2;
t156 = t207 ^ 2 + t212;
t210 = pkin(3) ^ 2;
t261 = t190 * t365;
t269 = 0.2e1 * pkin(2) * pkin(3);
t287 = t182 * t191;
t112 = t179 * t287 - t274;
t295 = t180 * t196;
t91 = t112 * t190 + t179 * t295;
t115 = t181 * t287 + t275;
t94 = t115 * t190 + t181 * t295;
t64 = (t204 * t91 + t234 * t94) * t82;
t353 = (-t207 * t261 + (t176 * t210 + t196 * t269 + t156) * t64) * t64;
t352 = t190 * pkin(2);
t351 = t192 * pkin(2);
t350 = t194 * pkin(2);
t260 = t192 * t364;
t285 = t182 * t193;
t113 = t179 * t285 - t272;
t293 = t180 * t198;
t92 = t113 * t192 + t179 * t293;
t116 = t181 * t285 + t273;
t95 = t116 * t192 + t181 * t293;
t65 = (t204 * t92 + t233 * t95) * t83;
t349 = (-t207 * t260 + (t177 * t210 + t198 * t269 + t156) * t65) * t65;
t259 = t194 * t363;
t283 = t182 * t195;
t114 = t179 * t283 - t270;
t291 = t180 * t200;
t93 = t114 * t194 + t179 * t291;
t117 = t181 * t283 + t271;
t96 = t117 * t194 + t181 * t291;
t66 = (t204 * t93 + t232 * t96) * t84;
t348 = (-t207 * t259 + (t178 * t210 + t200 * t269 + t156) * t66) * t66;
t347 = mrSges(3,1) * t190;
t346 = mrSges(3,1) * t192;
t345 = mrSges(3,1) * t194;
t341 = t164 * t94;
t340 = t165 * t95;
t339 = t166 * t96;
t338 = t167 * t94;
t337 = t168 * t95;
t336 = t169 * t96;
t335 = t197 * t64;
t334 = t199 * t65;
t333 = t201 * t66;
t332 = t207 * t64;
t331 = t207 * t65;
t330 = t207 * t66;
t329 = t211 * t76;
t328 = t211 * t77;
t327 = t211 * t78;
t326 = t211 * t79;
t325 = t211 * t80;
t324 = t211 * t81;
t246 = mrSges(3,2) * t196 + t347;
t97 = t182 * t249 - t246 * t300;
t323 = t211 * t97;
t245 = mrSges(3,2) * t198 + t346;
t98 = t182 * t248 - t245 * t298;
t322 = t211 * t98;
t244 = mrSges(3,2) * t200 + t345;
t99 = t182 * t247 - t244 * t296;
t321 = t211 * t99;
t170 = m(3) * pkin(2) + mrSges(2,1);
t124 = t249 + t170;
t157 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t100 = t124 * t197 + t157 * t191;
t320 = t100 * t180;
t125 = t248 + t170;
t101 = t125 * t199 + t157 * t193;
t319 = t101 * t180;
t126 = t247 + t170;
t102 = t126 * t201 + t157 * t195;
t318 = t102 * t180;
t118 = -t159 * t196 - t190 * t160;
t317 = t118 * t211;
t119 = -t159 * t198 - t192 * t160;
t316 = t119 * t211;
t120 = -t159 * t200 - t194 * t160;
t315 = t120 * t211;
t313 = t127 * t182;
t311 = t128 * t182;
t309 = t129 * t182;
t305 = mrSges(3,2) * t180 * t158;
t304 = t157 * t181;
t303 = t179 * t180;
t301 = t180 * t190;
t299 = t180 * t192;
t297 = t180 * t194;
t294 = t180 * t197;
t292 = t180 * t199;
t290 = t180 * t201;
t282 = t182 * t211;
t278 = t191 * t196;
t277 = t193 * t198;
t276 = t195 * t200;
t265 = -0.2e1 * t202;
t258 = t164 * t329;
t257 = t165 * t328;
t256 = t166 * t327;
t255 = t167 * t329;
t254 = t168 * t328;
t253 = t169 * t327;
t252 = t190 * t332;
t251 = t192 * t331;
t250 = t194 * t330;
t242 = t130 * t181 - t354;
t240 = t131 * t181 - t354;
t238 = t132 * t181 - t354;
t22 = t252 - t365;
t10 = (((t182 * t58 + t294 * t64) * t360 + ((-t261 + t332) * t191 + pkin(2) * t335) * t295 + t182 * t22) * t64 + (t58 * t294 + (t176 * t182 - t278 * t301 - t182) * t64) * t365) * t82;
t109 = pkin(3) * t278 + t133;
t13 = t82 * t282 * t353 + (-t182 * t252 + (-t109 * t301 + t182 * (pkin(2) * t196 + t360)) * t58) / (t109 * t295 + (pkin(2) + t357) * t288) * t58;
t142 = t157 * t354;
t16 = (-t196 * t353 - (pkin(2) * t58 - t196 * t22) * t365) * t82;
t216 = Ifges(3,1) + Ifges(2,3) + (pkin(6) ^ 2 + t212) * m(3) + 0.2e1 * pkin(6) * mrSges(3,3);
t88 = -t183 * t176 + 0.2e1 * (Ifges(3,4) * t190 + t203) * t196 + t190 * t265 + t216;
t4 = -t16 * t320 - t88 * t10 - t118 * t13 + ((t190 * t362 + t196 * t361) * t58 + (t203 * t190 + t367) * t64) * t378 + (-t124 * t379 - t130 * t304 + t142) * t197 + t191 * (t124 * t242 - t157 * t379);
t215 = t379 * t191 + t242 * t197;
t61 = t64 ^ 2;
t7 = -t97 * t16 - t118 * t10 - Ifges(3,3) * t13 + (pkin(2) * t347 + t367) * t61 + ((t180 * t243 - t313) * mrSges(3,1) + t215 * mrSges(3,2)) * t196 + t190 * (t305 + (t130 * t303 + t313) * mrSges(3,2) + t215 * mrSges(3,1));
t231 = t329 * t7 - t94 * t4;
t23 = t251 - t364;
t11 = (((t182 * t59 + t292 * t65) * t359 + ((-t260 + t331) * t193 + pkin(2) * t334) * t293 + t182 * t23) * t65 + (t59 * t292 + (t177 * t182 - t277 * t299 - t182) * t65) * t364) * t83;
t110 = pkin(3) * t277 + t134;
t14 = t83 * t282 * t349 + (-t182 * t251 + (-t110 * t299 + t182 * (pkin(2) * t198 + t359)) * t59) / (t110 * t293 + (pkin(2) + t356) * t286) * t59;
t17 = (-t198 * t349 - (pkin(2) * t59 - t198 * t23) * t364) * t83;
t89 = -t183 * t177 + 0.2e1 * (Ifges(3,4) * t192 + t203) * t198 + t192 * t265 + t216;
t5 = -t17 * t319 - t89 * t11 - t119 * t14 + ((t192 * t362 + t198 * t361) * t59 + (t203 * t192 + t368) * t65) * t377 + (-t125 * t380 - t131 * t304 + t142) * t199 + t193 * (t125 * t240 - t157 * t380);
t214 = t380 * t193 + t240 * t199;
t62 = t65 ^ 2;
t8 = -t98 * t17 - t119 * t11 - Ifges(3,3) * t14 + (pkin(2) * t346 + t368) * t62 + ((t180 * t241 - t311) * mrSges(3,1) + t214 * mrSges(3,2)) * t198 + t192 * (t305 + (t131 * t303 + t311) * mrSges(3,2) + t214 * mrSges(3,1));
t230 = t328 * t8 - t95 * t5;
t24 = t250 - t363;
t12 = (((t182 * t60 + t290 * t66) * t358 + ((-t259 + t330) * t195 + pkin(2) * t333) * t291 + t182 * t24) * t66 + (t60 * t290 + (t178 * t182 - t276 * t297 - t182) * t66) * t363) * t84;
t111 = pkin(3) * t276 + t135;
t15 = t84 * t282 * t348 + (-t182 * t250 + (-t111 * t297 + t182 * (pkin(2) * t200 + t358)) * t60) / (t111 * t291 + (pkin(2) + t355) * t284) * t60;
t18 = (-t200 * t348 - (pkin(2) * t60 - t200 * t24) * t363) * t84;
t90 = -t183 * t178 + 0.2e1 * (Ifges(3,4) * t194 + t203) * t200 + t194 * t265 + t216;
t6 = -t18 * t318 - t90 * t12 - t120 * t15 + ((t194 * t362 + t200 * t361) * t60 + (t203 * t194 + t369) * t66) * t376 + (-t126 * t381 - t132 * t304 + t142) * t201 + t195 * (t126 * t238 - t157 * t381);
t213 = t381 * t195 + t238 * t201;
t63 = t66 ^ 2;
t9 = -t99 * t18 - t120 * t12 - Ifges(3,3) * t15 + (pkin(2) * t345 + t369) * t63 + ((t180 * t239 - t309) * mrSges(3,1) + t213 * mrSges(3,2)) * t200 + t194 * (t305 + (t132 * t303 + t309) * mrSges(3,2) + t213 * mrSges(3,1));
t229 = t327 * t9 - t96 * t6;
t228 = Ifges(3,3) * t329 - t118 * t94;
t227 = Ifges(3,3) * t328 - t119 * t95;
t226 = Ifges(3,3) * t327 - t120 * t96;
t225 = t317 * t76 - t88 * t94;
t224 = t316 * t77 - t89 * t95;
t223 = t315 * t78 - t90 * t96;
t222 = pkin(3) * t301 - t133 * t182;
t221 = pkin(3) * t299 - t134 * t182;
t220 = pkin(3) * t297 - t135 * t182;
t219 = t320 * t94 - t323 * t76;
t218 = t319 * t95 - t322 * t77;
t217 = t318 * t96 - t321 * t78;
t189 = xDDP(1);
t188 = xDDP(2);
t187 = xDDP(3);
t174 = m(1) + m(2) + m(3);
t87 = t138 * t181 + t179 * t220;
t86 = t137 * t181 + t179 * t221;
t85 = t136 * t181 + t179 * t222;
t75 = -t117 * t358 - t138 * t179 * t200 + (pkin(2) * t297 + t200 * t220) * t181;
t74 = -t116 * t359 - t137 * t179 * t198 + (pkin(2) * t299 + t198 * t221) * t181;
t73 = -t115 * t360 - t136 * t179 * t196 + (pkin(2) * t301 + t196 * t222) * t181;
t72 = (t114 * t166 + t169 * t296) * t358 + (t105 * t169 - t166 * t87) * t200 + (-t166 * t303 + t169 * t182) * t350;
t71 = (t113 * t165 + t168 * t298) * t359 + (t104 * t168 - t165 * t86) * t198 + (-t165 * t303 + t168 * t182) * t351;
t70 = (t112 * t164 + t167 * t300) * t360 + (t103 * t167 - t164 * t85) * t196 + (-t164 * t303 + t167 * t182) * t352;
t69 = -(t114 * t169 - t166 * t296) * t358 + (t105 * t166 + t169 * t87) * t200 + (t166 * t182 + t169 * t303) * t350;
t68 = -(t113 * t168 - t165 * t298) * t359 + (t104 * t165 + t168 * t86) * t198 + (t165 * t182 + t168 * t303) * t351;
t67 = -(t112 * t167 - t164 * t300) * t360 + (t103 * t164 + t167 * t85) * t196 + (t164 * t182 + t167 * t303) * t352;
t57 = t60 ^ 2;
t56 = t59 ^ 2;
t55 = t58 ^ 2;
t54 = (Ifges(3,3) * t324 + t120 * t93 + t75 * t99) * t84;
t53 = (Ifges(3,3) * t325 + t119 * t92 + t74 * t98) * t83;
t52 = (Ifges(3,3) * t326 + t118 * t91 + t73 * t97) * t82;
t51 = (t174 * t75 + t318 * t93 + t321 * t81) * t84;
t50 = (t174 * t74 + t319 * t92 + t322 * t80) * t83;
t49 = (t174 * t73 + t320 * t91 + t323 * t79) * t82;
t48 = (t315 * t81 + t318 * t75 + t90 * t93) * t84;
t47 = (t316 * t80 + t319 * t74 + t89 * t92) * t83;
t46 = (t317 * t79 + t320 * t73 + t88 * t91) * t82;
t45 = (-t166 * t226 + t72 * t99) * t84;
t44 = (-t165 * t227 + t71 * t98) * t83;
t43 = (-t164 * t228 + t70 * t97) * t82;
t42 = (t169 * t226 + t69 * t99) * t84;
t41 = (t168 * t227 + t68 * t98) * t83;
t40 = (t167 * t228 + t67 * t97) * t82;
t39 = (t166 * t217 + t174 * t72) * t84;
t38 = (t165 * t218 + t174 * t71) * t83;
t37 = (t164 * t219 + t174 * t70) * t82;
t36 = (-t169 * t217 + t174 * t69) * t84;
t35 = (-t168 * t218 + t174 * t68) * t83;
t34 = (-t167 * t219 + t174 * t67) * t82;
t33 = (-t166 * t223 + t318 * t72) * t84;
t32 = (-t165 * t224 + t319 * t71) * t83;
t31 = (-t164 * t225 + t320 * t70) * t82;
t30 = (t169 * t223 + t318 * t69) * t84;
t29 = (t168 * t224 + t319 * t68) * t83;
t28 = (t167 * t225 + t320 * t67) * t82;
t3 = -t99 * t15 - t57 * t182 * t244 + (-t102 * t12 + (-t170 * t63 - t247 * (t63 + t57)) * t195 + (t157 * t66 + t244 * t376) * t333) * t180 + (-t18 - t129) * t174;
t2 = -t98 * t14 - t56 * t182 * t245 + (-t101 * t11 + (-t170 * t62 - t248 * (t62 + t56)) * t193 + (t157 * t65 + t245 * t377) * t334) * t180 + (-t17 - t128) * t174;
t1 = -t97 * t13 - t55 * t182 * t246 + (-t100 * t10 + (-t170 * t61 - t249 * (t61 + t55)) * t191 + (t157 * t64 + t246 * t378) * t335) * t180 + (-t16 - t127) * t174;
t19 = [(-g(1) + t189) * m(4) + ((-t256 * t42 + t30 * t339 + t36 * t72) * t188 + (t30 * t93 + t324 * t42 + t36 * t75) * t187 + (t36 * t189 + t3) * t69 + ((-t30 * t96 + t327 * t42) * t189 + t229) * t169) * t84 + ((-t257 * t41 + t29 * t340 + t35 * t71) * t188 + (t29 * t92 + t325 * t41 + t35 * t74) * t187 + (t35 * t189 + t2) * t68 + ((-t29 * t95 + t328 * t41) * t189 + t230) * t168) * t83 + ((-t258 * t40 + t28 * t341 + t34 * t70) * t188 + (t28 * t91 + t326 * t40 + t34 * t73) * t187 + (t34 * t189 + t1) * t67 + ((-t28 * t94 + t329 * t40) * t189 + t231) * t167) * t82; (-g(2) + t188) * m(4) + ((t253 * t45 - t33 * t336 + t39 * t69) * t189 + (t324 * t45 + t33 * t93 + t39 * t75) * t187 + (t39 * t188 + t3) * t72 + ((-t327 * t45 + t33 * t96) * t188 - t229) * t166) * t84 + ((t254 * t44 - t32 * t337 + t38 * t68) * t189 + (t32 * t92 + t325 * t44 + t38 * t74) * t187 + (t38 * t188 + t2) * t71 + ((t32 * t95 - t328 * t44) * t188 - t230) * t165) * t83 + ((t255 * t43 - t31 * t338 + t37 * t67) * t189 + (t31 * t91 + t326 * t43 + t37 * t73) * t187 + (t37 * t188 + t1) * t70 + ((t31 * t94 - t329 * t43) * t188 - t231) * t164) * t82; (-g(3) + t187) * m(4) + ((t253 * t54 - t336 * t48 + t51 * t69) * t189 + (-t256 * t54 + t339 * t48 + t51 * t72) * t188 + (t324 * t54 + t48 * t93 + t51 * t75) * t187 + t75 * t3 + t93 * t6 + t9 * t324) * t84 + ((t254 * t53 - t337 * t47 + t50 * t68) * t189 + (-t257 * t53 + t340 * t47 + t50 * t71) * t188 + (t325 * t53 + t47 * t92 + t50 * t74) * t187 + t74 * t2 + t92 * t5 + t8 * t325) * t83 + ((t255 * t52 - t338 * t46 + t49 * t67) * t189 + (-t258 * t52 + t341 * t46 + t49 * t70) * t188 + (t326 * t52 + t46 * t91 + t49 * t73) * t187 + t73 * t1 + t91 * t4 + t7 * t326) * t82;];
tauX  = t19;
