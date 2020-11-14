% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR8V2G2A0
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
% Datum: 2020-08-06 17:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:49:19
% EndTime: 2020-08-06 17:49:28
% DurationCPUTime: 8.66s
% Computational Cost: add. (39108->478), mult. (87021->871), div. (5112->7), fcn. (84471->22), ass. (0->316)
t194 = sin(qJ(2,1));
t200 = cos(qJ(2,1));
t206 = pkin(7) + pkin(6);
t135 = pkin(2) * t194 - t200 * t206;
t179 = sin(pkin(4));
t181 = cos(pkin(4));
t193 = sin(qJ(3,1));
t280 = t181 * t193;
t105 = pkin(3) * t280 + t135 * t179;
t199 = cos(qJ(3,1));
t295 = t179 * t194;
t177 = t199 ^ 2;
t353 = pkin(3) * t177;
t84 = 0.1e1 / (pkin(2) * t280 + t105 * t199 + t295 * t353);
t192 = sin(qJ(2,2));
t198 = cos(qJ(2,2));
t134 = pkin(2) * t192 - t198 * t206;
t191 = sin(qJ(3,2));
t282 = t181 * t191;
t104 = pkin(3) * t282 + t134 * t179;
t197 = cos(qJ(3,2));
t297 = t179 * t192;
t176 = t197 ^ 2;
t354 = pkin(3) * t176;
t83 = 0.1e1 / (pkin(2) * t282 + t104 * t197 + t297 * t354);
t190 = sin(qJ(2,3));
t196 = cos(qJ(2,3));
t133 = pkin(2) * t190 - t196 * t206;
t189 = sin(qJ(3,3));
t284 = t181 * t189;
t103 = pkin(3) * t284 + t133 * t179;
t195 = cos(qJ(3,3));
t299 = t179 * t190;
t175 = t195 ^ 2;
t355 = pkin(3) * t175;
t82 = 0.1e1 / (pkin(2) * t284 + t103 * t195 + t299 * t355);
t185 = legFrame(1,2);
t165 = sin(t185);
t168 = cos(t185);
t129 = t165 * g(1) + t168 * g(2);
t132 = t168 * g(1) - t165 * g(2);
t178 = sin(pkin(8));
t157 = g(3) * t178;
t180 = cos(pkin(8));
t237 = t132 * t180 - t157;
t379 = t129 * t179 + t237 * t181;
t184 = legFrame(2,2);
t164 = sin(t184);
t167 = cos(t184);
t128 = t164 * g(1) + t167 * g(2);
t131 = t167 * g(1) - t164 * g(2);
t239 = t131 * t180 - t157;
t378 = t128 * t179 + t239 * t181;
t183 = legFrame(3,2);
t163 = sin(t183);
t166 = cos(t183);
t127 = t163 * g(1) + t166 * g(2);
t130 = t166 * g(1) - t163 * g(2);
t241 = t130 * t180 - t157;
t377 = t127 * t179 + t241 * t181;
t203 = xDP(3);
t210 = 0.1e1 / pkin(3);
t204 = xDP(2);
t205 = xDP(1);
t233 = t163 * t204 - t166 * t205;
t136 = pkin(2) * t196 + t190 * t206;
t278 = t181 * t196;
t288 = t180 * t181;
t350 = t195 * pkin(3);
t76 = (t178 * t190 - t180 * t278) * t350 - t136 * t288 + t178 * t133;
t302 = t178 * t181;
t79 = (t178 * t278 + t180 * t190) * t350 + t136 * t302 + t133 * t180;
t58 = (t203 * t76 + t233 * t79) * t82 * t210;
t376 = -0.2e1 * t58;
t232 = t164 * t204 - t167 * t205;
t137 = pkin(2) * t198 + t192 * t206;
t277 = t181 * t198;
t349 = t197 * pkin(3);
t77 = (t178 * t192 - t180 * t277) * t349 - t137 * t288 + t178 * t134;
t80 = (t178 * t277 + t180 * t192) * t349 + t137 * t302 + t134 * t180;
t59 = (t203 * t77 + t232 * t80) * t83 * t210;
t375 = -0.2e1 * t59;
t231 = t165 * t204 - t168 * t205;
t138 = pkin(2) * t200 + t194 * t206;
t276 = t181 * t200;
t348 = t199 * pkin(3);
t78 = (t178 * t194 - t180 * t276) * t348 - t138 * t288 + t178 * t135;
t81 = (t178 * t276 + t180 * t194) * t348 + t138 * t302 + t135 * t180;
t60 = (t203 * t78 + t231 * t81) * t84 * t210;
t374 = -0.2e1 * t60;
t248 = t195 * mrSges(3,1) - mrSges(3,2) * t189;
t247 = t197 * mrSges(3,1) - mrSges(3,2) * t191;
t246 = t199 * mrSges(3,1) - mrSges(3,2) * t193;
t182 = Ifges(3,1) - Ifges(3,2);
t201 = pkin(2) * mrSges(3,2);
t364 = -2 * Ifges(3,4);
t367 = Ifges(3,4) + t177 * t364 + (-t182 * t193 + t201) * t199;
t366 = Ifges(3,4) + t176 * t364 + (-t182 * t191 + t201) * t197;
t365 = Ifges(3,4) + t175 * t364 + (-t182 * t189 + t201) * t195;
t202 = mrSges(3,1) * pkin(2);
t363 = pkin(3) * t58;
t362 = pkin(3) * t59;
t361 = pkin(3) * t60;
t158 = mrSges(3,2) * pkin(6) - Ifges(3,6);
t360 = -t158 / 0.2e1;
t159 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t359 = t159 / 0.2e1;
t358 = pkin(2) * t189;
t357 = pkin(2) * t191;
t356 = pkin(2) * t193;
t352 = g(3) * t180;
t211 = pkin(2) ^ 2;
t155 = t206 ^ 2 + t211;
t209 = pkin(3) ^ 2;
t260 = t189 * t363;
t268 = 0.2e1 * pkin(2) * pkin(3);
t283 = t181 * t190;
t112 = t178 * t283 - t180 * t196;
t294 = t179 * t195;
t91 = t189 * t112 + t178 * t294;
t115 = t178 * t196 + t180 * t283;
t287 = t180 * t195;
t94 = -t189 * t115 - t179 * t287;
t64 = (t203 * t94 + t233 * t91) * t82;
t351 = (-t206 * t260 + (t175 * t209 + t195 * t268 + t155) * t64) * t64;
t259 = t191 * t362;
t281 = t181 * t192;
t113 = t178 * t281 - t180 * t198;
t292 = t179 * t197;
t92 = t191 * t113 + t178 * t292;
t116 = t178 * t198 + t180 * t281;
t286 = t180 * t197;
t95 = -t191 * t116 - t179 * t286;
t65 = (t203 * t95 + t232 * t92) * t83;
t347 = (-t206 * t259 + (t176 * t209 + t197 * t268 + t155) * t65) * t65;
t258 = t193 * t361;
t279 = t181 * t194;
t114 = t178 * t279 - t180 * t200;
t290 = t179 * t199;
t93 = t193 * t114 + t178 * t290;
t117 = t178 * t200 + t180 * t279;
t285 = t180 * t199;
t96 = -t193 * t117 - t179 * t285;
t66 = (t203 * t96 + t231 * t93) * t84;
t346 = (-t206 * t258 + (t177 * t209 + t199 * t268 + t155) * t66) * t66;
t342 = t163 * t91;
t341 = t164 * t92;
t340 = t165 * t93;
t339 = t166 * t91;
t338 = t167 * t92;
t337 = t168 * t93;
t336 = t189 * mrSges(3,1);
t335 = t191 * mrSges(3,1);
t334 = t193 * mrSges(3,1);
t333 = t196 * t64;
t332 = t198 * t65;
t331 = t200 * t66;
t330 = t210 * t76;
t329 = t210 * t77;
t328 = t210 * t78;
t327 = t210 * t79;
t326 = t210 * t80;
t325 = t210 * t81;
t245 = t195 * mrSges(3,2) + t336;
t97 = t248 * t181 - t245 * t299;
t324 = t210 * t97;
t244 = t197 * mrSges(3,2) + t335;
t98 = t247 * t181 - t244 * t297;
t323 = t210 * t98;
t243 = t199 * mrSges(3,2) + t334;
t99 = t246 * t181 - t243 * t295;
t322 = t210 * t99;
t321 = t64 * t206;
t320 = t65 * t206;
t319 = t66 * t206;
t169 = m(3) * pkin(2) + mrSges(2,1);
t124 = t248 + t169;
t156 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t100 = t124 * t196 + t190 * t156;
t318 = t100 * t179;
t125 = t247 + t169;
t101 = t125 * t198 + t192 * t156;
t317 = t101 * t179;
t126 = t246 + t169;
t102 = t126 * t200 + t194 * t156;
t316 = t102 * t179;
t118 = -t158 * t195 - t189 * t159;
t315 = t118 * t210;
t119 = -t158 * t197 - t191 * t159;
t314 = t119 * t210;
t120 = -t158 * t199 - t193 * t159;
t313 = t120 * t210;
t311 = t127 * t181;
t309 = t128 * t181;
t307 = t129 * t181;
t303 = mrSges(3,2) * t157 * t179;
t301 = t179 * t180;
t300 = t179 * t189;
t298 = t179 * t191;
t296 = t179 * t193;
t293 = t179 * t196;
t291 = t179 * t198;
t289 = t179 * t200;
t275 = t181 * t210;
t271 = t190 * t195;
t270 = t192 * t197;
t269 = t194 * t199;
t264 = -0.2e1 * t201;
t257 = t163 * t327;
t256 = t164 * t326;
t255 = t165 * t325;
t254 = t166 * t327;
t253 = t167 * t326;
t252 = t168 * t325;
t251 = t189 * t321;
t250 = t191 * t320;
t249 = t193 * t319;
t242 = t178 * t130 + t352;
t240 = t178 * t131 + t352;
t238 = t178 * t132 + t352;
t22 = t251 - t363;
t10 = (((t181 * t58 + t64 * t293) * t355 + ((-t260 + t321) * t190 + pkin(2) * t333) * t294 + t181 * t22) * t64 + (t58 * t293 + (t175 * t181 - t271 * t300 - t181) * t64) * t363) * t82;
t109 = pkin(3) * t271 + t133;
t13 = t82 * t275 * t351 + (-t181 * t251 + (-t109 * t300 + t181 * (pkin(2) * t195 + t355)) * t58) / (t109 * t294 + (pkin(2) + t350) * t284) * t58;
t16 = (-t195 * t351 - (pkin(2) * t58 - t22 * t195) * t363) * t82;
t215 = Ifges(3,1) + Ifges(2,3) + (pkin(6) ^ 2 + t211) * m(3) + 0.2e1 * pkin(6) * mrSges(3,3);
t88 = -t182 * t175 + 0.2e1 * (Ifges(3,4) * t189 + t202) * t195 + t189 * t264 + t215;
t4 = -t16 * t318 - t88 * t10 - t118 * t13 + ((t360 * t189 + t359 * t195) * t58 + (t202 * t189 + t365) * t64) * t376 + (-t124 * t377 - t242 * t156) * t196 + t190 * (t242 * t124 - t156 * t377);
t214 = t377 * t190 + t242 * t196;
t61 = t64 ^ 2;
t7 = -t97 * t16 - t118 * t10 - Ifges(3,3) * t13 + (pkin(2) * t336 + t365) * t61 + ((t241 * t179 - t311) * mrSges(3,1) + t214 * mrSges(3,2)) * t195 + t189 * (t303 + (-t130 * t301 + t311) * mrSges(3,2) + t214 * mrSges(3,1));
t230 = t7 * t327 + t91 * t4;
t23 = t250 - t362;
t11 = (((t181 * t59 + t65 * t291) * t354 + ((-t259 + t320) * t192 + pkin(2) * t332) * t292 + t181 * t23) * t65 + (t59 * t291 + (t176 * t181 - t270 * t298 - t181) * t65) * t362) * t83;
t110 = pkin(3) * t270 + t134;
t14 = t83 * t275 * t347 + (-t181 * t250 + (-t110 * t298 + t181 * (pkin(2) * t197 + t354)) * t59) / (t110 * t292 + (pkin(2) + t349) * t282) * t59;
t17 = (-t197 * t347 - (pkin(2) * t59 - t23 * t197) * t362) * t83;
t89 = -t182 * t176 + 0.2e1 * (Ifges(3,4) * t191 + t202) * t197 + t191 * t264 + t215;
t5 = -t17 * t317 - t89 * t11 - t119 * t14 + ((t360 * t191 + t359 * t197) * t59 + (t202 * t191 + t366) * t65) * t375 + (-t125 * t378 - t240 * t156) * t198 + t192 * (t240 * t125 - t156 * t378);
t213 = t378 * t192 + t240 * t198;
t62 = t65 ^ 2;
t8 = -t98 * t17 - t119 * t11 - Ifges(3,3) * t14 + (pkin(2) * t335 + t366) * t62 + ((t239 * t179 - t309) * mrSges(3,1) + t213 * mrSges(3,2)) * t197 + t191 * (t303 + (-t131 * t301 + t309) * mrSges(3,2) + t213 * mrSges(3,1));
t229 = t8 * t326 + t92 * t5;
t24 = t249 - t361;
t12 = (((t181 * t60 + t66 * t289) * t353 + ((-t258 + t319) * t194 + pkin(2) * t331) * t290 + t181 * t24) * t66 + (t60 * t289 + (t177 * t181 - t269 * t296 - t181) * t66) * t361) * t84;
t111 = pkin(3) * t269 + t135;
t15 = t84 * t275 * t346 + (-t181 * t249 + (-t111 * t296 + t181 * (pkin(2) * t199 + t353)) * t60) / (t111 * t290 + (pkin(2) + t348) * t280) * t60;
t18 = (-t199 * t346 - (pkin(2) * t60 - t24 * t199) * t361) * t84;
t90 = -t182 * t177 + 0.2e1 * (Ifges(3,4) * t193 + t202) * t199 + t193 * t264 + t215;
t6 = -t18 * t316 - t90 * t12 - t120 * t15 + ((t360 * t193 + t359 * t199) * t60 + (t202 * t193 + t367) * t66) * t374 + (-t126 * t379 - t238 * t156) * t200 + t194 * (t238 * t126 - t156 * t379);
t212 = t379 * t194 + t238 * t200;
t63 = t66 ^ 2;
t9 = -t99 * t18 - t120 * t12 - Ifges(3,3) * t15 + (pkin(2) * t334 + t367) * t63 + ((t237 * t179 - t307) * mrSges(3,1) + t212 * mrSges(3,2)) * t199 + t193 * (t303 + (-t132 * t301 + t307) * mrSges(3,2) + t212 * mrSges(3,1));
t228 = t9 * t325 + t93 * t6;
t227 = Ifges(3,3) * t327 + t118 * t91;
t226 = Ifges(3,3) * t326 + t119 * t92;
t225 = Ifges(3,3) * t325 + t120 * t93;
t224 = t79 * t315 + t88 * t91;
t223 = t80 * t314 + t89 * t92;
t222 = t81 * t313 + t90 * t93;
t221 = pkin(3) * t300 - t133 * t181;
t220 = pkin(3) * t298 - t134 * t181;
t219 = pkin(3) * t296 - t135 * t181;
t218 = t91 * t318 + t79 * t324;
t217 = t92 * t317 + t80 * t323;
t216 = t93 * t316 + t81 * t322;
t188 = xDDP(1);
t187 = xDDP(2);
t186 = xDDP(3);
t173 = m(1) + m(2) + m(3);
t87 = -t178 * t138 + t219 * t180;
t86 = -t178 * t137 + t220 * t180;
t85 = -t178 * t136 + t221 * t180;
t75 = -t114 * t353 + t138 * t285 + (pkin(2) * t296 + t219 * t199) * t178;
t74 = -t113 * t354 + t137 * t286 + (pkin(2) * t298 + t220 * t197) * t178;
t73 = -t112 * t355 + t136 * t287 + (pkin(2) * t300 + t221 * t195) * t178;
t72 = (t117 * t168 + t165 * t295) * t353 + (t105 * t165 - t87 * t168) * t199 + (t181 * t165 - t168 * t301) * t356;
t71 = -(t117 * t165 - t168 * t295) * t353 + (t168 * t105 + t87 * t165) * t199 + (t165 * t301 + t168 * t181) * t356;
t70 = (t116 * t167 + t164 * t297) * t354 + (t104 * t164 - t86 * t167) * t197 + (t181 * t164 - t167 * t301) * t357;
t69 = -(t116 * t164 - t167 * t297) * t354 + (t167 * t104 + t86 * t164) * t197 + (t164 * t301 + t167 * t181) * t357;
t68 = (t115 * t166 + t163 * t299) * t355 + (t103 * t163 - t85 * t166) * t195 + (t181 * t163 - t166 * t301) * t358;
t67 = -(t115 * t163 - t166 * t299) * t355 + (t166 * t103 + t85 * t163) * t195 + (t163 * t301 + t166 * t181) * t358;
t57 = t60 ^ 2;
t56 = t59 ^ 2;
t55 = t58 ^ 2;
t54 = (Ifges(3,3) * t328 + t120 * t96 + t75 * t99) * t84;
t53 = (Ifges(3,3) * t329 + t119 * t95 + t74 * t98) * t83;
t52 = (Ifges(3,3) * t330 + t118 * t94 + t73 * t97) * t82;
t51 = (t173 * t75 + t96 * t316 + t78 * t322) * t84;
t50 = (t173 * t74 + t95 * t317 + t77 * t323) * t83;
t49 = (t173 * t73 + t94 * t318 + t76 * t324) * t82;
t48 = (t78 * t313 + t75 * t316 + t90 * t96) * t84;
t47 = (t77 * t314 + t74 * t317 + t89 * t95) * t83;
t46 = (t76 * t315 + t73 * t318 + t88 * t94) * t82;
t45 = (-t225 * t168 + t72 * t99) * t84;
t44 = (t225 * t165 + t71 * t99) * t84;
t43 = (-t226 * t167 + t70 * t98) * t83;
t42 = (t226 * t164 + t69 * t98) * t83;
t41 = (-t227 * t166 + t68 * t97) * t82;
t40 = (t227 * t163 + t67 * t97) * t82;
t39 = (-t216 * t168 + t173 * t72) * t84;
t38 = (t216 * t165 + t173 * t71) * t84;
t37 = (-t217 * t167 + t173 * t70) * t83;
t36 = (t217 * t164 + t173 * t69) * t83;
t35 = (-t218 * t166 + t173 * t68) * t82;
t34 = (t218 * t163 + t173 * t67) * t82;
t33 = (-t222 * t168 + t72 * t316) * t84;
t32 = (t222 * t165 + t71 * t316) * t84;
t31 = (-t223 * t167 + t70 * t317) * t83;
t30 = (t223 * t164 + t69 * t317) * t83;
t29 = (-t224 * t166 + t68 * t318) * t82;
t28 = (t224 * t163 + t67 * t318) * t82;
t3 = -t99 * t15 - t57 * t181 * t243 + (-t102 * t12 + (-t63 * t169 - t246 * (t63 + t57)) * t194 + (t66 * t156 + t243 * t374) * t331) * t179 + (-t18 - t129) * t173;
t2 = -t98 * t14 - t56 * t181 * t244 + (-t101 * t11 + (-t62 * t169 - t247 * (t62 + t56)) * t192 + (t65 * t156 + t244 * t375) * t332) * t179 + (-t17 - t128) * t173;
t1 = -t97 * t13 - t55 * t181 * t245 + (-t100 * t10 + (-t61 * t169 - t248 * (t61 + t55)) * t190 + (t64 * t156 + t245 * t376) * t333) * t179 + (-t16 - t127) * t173;
t19 = [(-g(1) + t188) * m(4) + ((t45 * t255 + t33 * t340 + t39 * t71) * t187 + (t45 * t328 + t33 * t96 + t39 * t75) * t186 + (t39 * t188 + t3) * t72 + ((-t45 * t325 - t33 * t93) * t188 - t228) * t168) * t84 + ((t43 * t256 + t31 * t341 + t37 * t69) * t187 + (t31 * t95 + t43 * t329 + t37 * t74) * t186 + (t37 * t188 + t2) * t70 + ((-t31 * t92 - t43 * t326) * t188 - t229) * t167) * t83 + ((t41 * t257 + t29 * t342 + t35 * t67) * t187 + (t29 * t94 + t41 * t330 + t35 * t73) * t186 + (t35 * t188 + t1) * t68 + ((-t29 * t91 - t41 * t327) * t188 - t230) * t166) * t82; (-g(2) + t187) * m(4) + ((-t44 * t252 - t32 * t337 + t38 * t72) * t188 + (t32 * t96 + t44 * t328 + t38 * t75) * t186 + (t38 * t187 + t3) * t71 + ((t32 * t93 + t44 * t325) * t187 + t228) * t165) * t84 + ((-t42 * t253 - t30 * t338 + t36 * t70) * t188 + (t30 * t95 + t42 * t329 + t36 * t74) * t186 + (t36 * t187 + t2) * t69 + ((t30 * t92 + t42 * t326) * t187 + t229) * t164) * t83 + ((-t40 * t254 - t28 * t339 + t34 * t68) * t188 + (t28 * t94 + t40 * t330 + t34 * t73) * t186 + (t34 * t187 + t1) * t67 + ((t28 * t91 + t40 * t327) * t187 + t230) * t163) * t82; (-g(3) + t186) * m(4) + ((-t54 * t252 - t48 * t337 + t51 * t72) * t188 + (t54 * t255 + t48 * t340 + t51 * t71) * t187 + (t54 * t328 + t48 * t96 + t51 * t75) * t186 + t75 * t3 + t96 * t6 + t9 * t328) * t84 + ((-t53 * t253 - t47 * t338 + t50 * t70) * t188 + (t53 * t256 + t47 * t341 + t50 * t69) * t187 + (t53 * t329 + t47 * t95 + t50 * t74) * t186 + t74 * t2 + t95 * t5 + t8 * t329) * t83 + ((-t52 * t254 - t46 * t339 + t49 * t68) * t188 + (t52 * t257 + t46 * t342 + t49 * t67) * t187 + (t52 * t330 + t46 * t94 + t49 * t73) * t186 + t73 * t1 + t94 * t4 + t7 * t330) * t82;];
tauX  = t19;
