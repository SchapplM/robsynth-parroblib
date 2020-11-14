% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRRR9V1G2A0
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
% Datum: 2020-08-06 18:53
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:51:31
% EndTime: 2020-08-06 18:51:38
% DurationCPUTime: 6.50s
% Computational Cost: add. (24843->621), mult. (32331->1033), div. (4665->10), fcn. (26568->36), ass. (0->373)
t202 = cos(pkin(7));
t160 = t202 * pkin(2);
t412 = -0.2e1 * t160;
t410 = 2 * pkin(1);
t411 = t410 / 0.2e1;
t231 = (m(2) + m(3));
t184 = pkin(1) * t231;
t147 = t160 + pkin(1);
t192 = pkin(7) + qJ(3,1);
t159 = cos(t192);
t395 = pkin(3) * t159;
t114 = t147 + t395;
t191 = pkin(7) + qJ(3,2);
t158 = cos(t191);
t396 = pkin(3) * t158;
t113 = t147 + t396;
t190 = pkin(7) + qJ(3,3);
t157 = cos(t190);
t397 = pkin(3) * t157;
t112 = t147 + t397;
t203 = (pkin(5) + qJ(2,3));
t140 = mrSges(3,1) * t203 - Ifges(3,5);
t204 = (pkin(5) + qJ(2,2));
t141 = mrSges(3,1) * t204 - Ifges(3,5);
t205 = (pkin(5) + qJ(2,1));
t142 = mrSges(3,1) * t205 - Ifges(3,5);
t137 = mrSges(3,2) * t203 - Ifges(3,6);
t138 = mrSges(3,2) * t204 - Ifges(3,6);
t139 = mrSges(3,2) * t205 - Ifges(3,6);
t409 = 0.2e1 * pkin(3);
t189 = t202 ^ 2;
t408 = 0.2e1 * t189;
t201 = sin(pkin(7));
t407 = -0.2e1 * t201;
t406 = 0.4e1 * t202;
t207 = (mrSges(2,3) + mrSges(3,3));
t405 = 2 * t207;
t404 = -4 * pkin(5) - 4 * pkin(6);
t403 = pkin(1) * mrSges(3,2);
t151 = 0.1e1 / t157;
t186 = pkin(6) + t203;
t168 = 1 / t186;
t221 = cos(qJ(1,3));
t227 = xDP(3);
t228 = xDP(2);
t229 = xDP(1);
t154 = sin(t190);
t208 = legFrame(3,2);
t174 = cos(t208);
t171 = sin(t208);
t215 = sin(qJ(1,3));
t337 = t171 * t215;
t97 = t154 * t174 - t157 * t337;
t331 = t174 * t215;
t98 = t154 * t171 + t157 * t331;
t67 = (t221 * t227 + (t228 * t97 + t229 * t98) * t151) * t168;
t65 = pkin(1) * t67;
t155 = sin(t191);
t209 = legFrame(2,2);
t172 = sin(t209);
t175 = cos(t209);
t217 = sin(qJ(1,2));
t329 = t175 * t217;
t100 = t155 * t172 + t158 * t329;
t152 = 0.1e1 / t158;
t187 = pkin(6) + t204;
t169 = 1 / t187;
t223 = cos(qJ(1,2));
t335 = t172 * t217;
t99 = t155 * t175 - t158 * t335;
t68 = (t223 * t227 + (t100 * t229 + t228 * t99) * t152) * t169;
t66 = pkin(1) * t68;
t402 = pkin(2) * mrSges(3,1);
t401 = pkin(2) * mrSges(3,2);
t156 = sin(t192);
t210 = legFrame(1,2);
t176 = cos(t210);
t173 = sin(t210);
t219 = sin(qJ(1,1));
t333 = t173 * t219;
t101 = t156 * t176 - t159 * t333;
t327 = t176 * t219;
t102 = t156 * t173 + t159 * t327;
t153 = 0.1e1 / t159;
t188 = pkin(6) + t205;
t170 = 1 / t188;
t225 = cos(qJ(1,1));
t69 = (t225 * t227 + (t101 * t228 + t102 * t229) * t153) * t170;
t64 = t69 * pkin(1);
t224 = cos(qJ(3,1));
t218 = sin(qJ(3,1));
t323 = t201 * t218;
t108 = 0.1e1 / (t202 * t224 - t323);
t345 = t108 * t170;
t244 = -pkin(3) / 0.2e1;
t389 = t224 * pkin(2);
t199 = t224 ^ 2;
t392 = t199 * pkin(3);
t123 = t392 + t389 / 0.2e1 + t244;
t274 = pkin(1) * t219 - t188 * t225;
t277 = t219 * t323;
t258 = pkin(2) * t277 + (t277 * t409 - t274) * t224;
t245 = pkin(2) / 0.2e1;
t338 = (pkin(3) * t224 + t245) * t218;
t161 = pkin(1) * t201;
t341 = (-pkin(3) * t218 + t161) * t224;
t84 = t274 * t323 + (t199 - 0.1e1) * t219 * pkin(3);
t96 = pkin(1) * t218 + (-pkin(3) + t389 + 0.2e1 * t392) * t201;
t50 = (-t123 * t333 + t176 * t338) * t408 + (t258 * t173 + t176 * t96) * t202 + t84 * t173 + t176 * t341;
t44 = t50 * t228 * t345;
t51 = (t123 * t327 + t173 * t338) * t408 + (t173 * t96 - t258 * t176) * t202 - t84 * t176 + t173 * t341;
t45 = t51 * t229 * t345;
t93 = t114 * t225 + t188 * t219;
t81 = t93 * t170 * t227;
t19 = t64 - 0.2e1 * t45 - 0.2e1 * t44 - 0.2e1 * t81;
t400 = mrSges(3,2) * t19;
t220 = cos(qJ(3,3));
t214 = sin(qJ(3,3));
t325 = t201 * t214;
t106 = 0.1e1 / (t202 * t220 - t325);
t349 = t106 * t168;
t391 = t220 * pkin(2);
t197 = t220 ^ 2;
t394 = t197 * pkin(3);
t121 = t394 + t391 / 0.2e1 + t244;
t276 = pkin(1) * t215 - t186 * t221;
t279 = t215 * t325;
t260 = pkin(2) * t279 + (t279 * t409 - t276) * t220;
t340 = (pkin(3) * t220 + t245) * t214;
t343 = (-pkin(3) * t214 + t161) * t220;
t82 = t276 * t325 + (t197 - 0.1e1) * t215 * pkin(3);
t94 = pkin(1) * t214 + (-pkin(3) + t391 + 0.2e1 * t394) * t201;
t46 = (-t121 * t337 + t174 * t340) * t408 + (t260 * t171 + t174 * t94) * t202 + t82 * t171 + t174 * t343;
t40 = t46 * t228 * t349;
t47 = (t121 * t331 + t171 * t340) * t408 + (t171 * t94 - t260 * t174) * t202 - t82 * t174 + t171 * t343;
t41 = t47 * t229 * t349;
t91 = t112 * t221 + t186 * t215;
t79 = t91 * t168 * t227;
t20 = t65 - 0.2e1 * t41 - 0.2e1 * t40 - 0.2e1 * t79;
t399 = mrSges(3,2) * t20;
t222 = cos(qJ(3,2));
t216 = sin(qJ(3,2));
t324 = t201 * t216;
t107 = 0.1e1 / (t202 * t222 - t324);
t347 = t107 * t169;
t390 = t222 * pkin(2);
t198 = t222 ^ 2;
t393 = t198 * pkin(3);
t122 = t393 + t390 / 0.2e1 + t244;
t275 = pkin(1) * t217 - t187 * t223;
t278 = t217 * t324;
t259 = pkin(2) * t278 + (t278 * t409 - t275) * t222;
t339 = (pkin(3) * t222 + t245) * t216;
t342 = (-pkin(3) * t216 + t161) * t222;
t83 = t275 * t324 + (t198 - 0.1e1) * t217 * pkin(3);
t95 = pkin(1) * t216 + (-pkin(3) + t390 + 0.2e1 * t393) * t201;
t48 = (-t122 * t335 + t175 * t339) * t408 + (t259 * t172 + t175 * t95) * t202 + t83 * t172 + t175 * t342;
t42 = t48 * t228 * t347;
t49 = (t122 * t329 + t172 * t339) * t408 + (t172 * t95 - t259 * t175) * t202 - t83 * t175 + t172 * t342;
t43 = t49 * t229 * t347;
t92 = t113 * t223 + t187 * t217;
t80 = t92 * t169 * t227;
t21 = t66 - 0.2e1 * t43 - 0.2e1 * t42 - 0.2e1 * t80;
t398 = mrSges(3,2) * t21;
t206 = Ifges(3,1) - Ifges(3,2);
t388 = 2 * Ifges(2,4) - 2 * Ifges(3,4);
t387 = mrSges(3,1) * t201;
t386 = mrSges(3,2) * t214;
t385 = mrSges(3,2) * t216;
t384 = mrSges(3,2) * t218;
t383 = Ifges(3,4) * t214;
t382 = Ifges(3,4) * t216;
t381 = Ifges(3,4) * t218;
t380 = t106 * t46;
t379 = t106 * t47;
t378 = t107 * t48;
t377 = t107 * t49;
t376 = t108 * t50;
t375 = t108 * t51;
t250 = 0.1e1 / pkin(3);
t88 = (t171 * t229 + t174 * t228) * t250 * t151;
t374 = t88 ^ 2 * t151;
t373 = t151 * t97;
t372 = t151 * t98;
t89 = (t172 * t229 + t175 * t228) * t250 * t152;
t371 = t89 ^ 2 * t152;
t370 = t152 * t99;
t90 = (t173 * t229 + t176 * t228) * t250 * t153;
t369 = t90 ^ 2 * t153;
t368 = t168 * t97;
t367 = t168 * t98;
t366 = t169 * t99;
t365 = t189 * t67;
t364 = t189 * t68;
t363 = t189 * t69;
t362 = t197 * Ifges(3,4);
t361 = t198 * Ifges(3,4);
t360 = t199 * Ifges(3,4);
t359 = t201 * t67;
t358 = t201 * t68;
t357 = t201 * t69;
t162 = t214 * mrSges(3,1);
t163 = t216 * mrSges(3,1);
t164 = t218 * mrSges(3,1);
t356 = m(3) * pkin(2) + mrSges(2,1);
t355 = t100 * t152;
t354 = t100 * t169;
t353 = t101 * t153;
t352 = t101 * t170;
t351 = t102 * t153;
t350 = t102 * t170;
t348 = t106 * t231;
t346 = t107 * t231;
t344 = t108 * t231;
t336 = t171 * t250;
t334 = t172 * t250;
t332 = t173 * t250;
t330 = t174 * t250;
t328 = t175 * t250;
t326 = t176 * t250;
t322 = t206 * t197;
t321 = t206 * t198;
t320 = t206 * t199;
t319 = t206 * t214;
t318 = t206 * t216;
t317 = t206 * t218;
t316 = mrSges(3,2) * t220 + t162;
t315 = mrSges(3,2) * t222 + t163;
t314 = mrSges(3,2) * t224 + t164;
t313 = -2 * t403;
t312 = -0.2e1 * t161;
t310 = mrSges(3,1) * t65;
t309 = mrSges(3,1) * t66;
t308 = mrSges(3,1) * t64;
t307 = 0.4e1 * t383;
t306 = 0.4e1 * t382;
t305 = 0.4e1 * t381;
t303 = mrSges(3,2) * t161;
t302 = pkin(2) * t386;
t301 = pkin(2) * t385;
t300 = pkin(2) * t384;
t299 = -mrSges(1,2) + t207;
t28 = t41 + t40 + t79;
t29 = t43 + t42 + t80;
t30 = t45 + t44 + t81;
t298 = t67 * t362;
t297 = t68 * t361;
t296 = t69 * t360;
t295 = -t403 / 0.2e1;
t294 = -t402 / 0.2e1;
t293 = -t402 / 0.4e1;
t292 = m(3) * pkin(5) + t207;
t269 = mrSges(3,1) * t220 - t386;
t263 = (-t269 - t356) * t202 + (mrSges(2,2) + t316) * t201;
t76 = -t184 + t263;
t291 = t76 * t349;
t268 = mrSges(3,1) * t222 - t385;
t262 = (-t268 - t356) * t202 + (mrSges(2,2) + t315) * t201;
t77 = -t184 + t262;
t290 = t77 * t347;
t267 = mrSges(3,1) * t224 - t384;
t261 = (-t267 - t356) * t202 + (mrSges(2,2) + t314) * t201;
t78 = -t184 + t261;
t289 = t78 * t345;
t288 = t154 * t374;
t287 = t155 * t371;
t286 = t156 * t369;
t73 = (-t137 * t220 - t140 * t214) * t202 + t201 * (t137 * t214 - t140 * t220);
t285 = t221 * t250 * t73;
t74 = (-t138 * t222 - t141 * t216) * t202 + t201 * (t138 * t216 - t141 * t222);
t284 = t223 * t250 * t74;
t75 = (-t139 * t224 - t142 * t218) * t202 + t201 * (t139 * t218 - t142 * t224);
t283 = t225 * t250 * t75;
t282 = (t302 - t206) * t359;
t281 = (t301 - t206) * t358;
t280 = (t300 - t206) * t357;
t273 = t322 * t359;
t272 = t321 * t358;
t271 = t320 * t357;
t251 = pkin(2) ^ 2;
t270 = t251 * m(3) - Ifges(2,1) + Ifges(2,2) + t206;
t118 = g(1) * t174 - g(2) * t171;
t266 = g(3) * t221 + t118 * t215;
t119 = g(1) * t175 - g(2) * t172;
t265 = g(3) * t223 + t119 * t217;
t120 = g(1) * t176 - g(2) * t173;
t264 = g(3) * t225 + t120 * t219;
t257 = m(2) * qJ(2,1) + m(3) * t205 + t299;
t256 = m(2) * qJ(2,2) + m(3) * t204 + t299;
t255 = m(2) * qJ(2,3) + m(3) * t203 + t299;
t252 = pkin(1) ^ 2;
t254 = (2 * mrSges(3,3) * pkin(5)) + (t231 * t252) + Ifges(2,1) + Ifges(3,2) + Ifges(1,3);
t241 = 0.2e1 * pkin(7);
t249 = pkin(3) ^ 2;
t253 = -t251 * cos(t241) - (2 * pkin(6) ^ 2) - t249 - t251 - (2 * t252) + ((-4 * pkin(6) - 2 * pkin(5)) * pkin(5));
t248 = qJ(2,1) ^ 2;
t247 = qJ(2,2) ^ 2;
t246 = qJ(2,3) ^ 2;
t243 = -pkin(5) / 0.2e1;
t242 = -pkin(5) / 0.4e1;
t234 = -Ifges(3,4) / 0.2e1;
t233 = Ifges(3,5) / 0.4e1;
t232 = Ifges(3,6) / 0.2e1;
t213 = xDDP(1);
t212 = xDDP(2);
t211 = xDDP(3);
t196 = mrSges(3,1) * t410;
t194 = 0.2e1 * t402;
t193 = -0.2e1 * t401;
t183 = t402 / 0.4e1;
t182 = -t401 / 0.2e1;
t181 = -t401 / 0.4e1;
t177 = Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1;
t150 = 0.2e1 * t192;
t149 = 0.2e1 * t191;
t148 = 0.2e1 * t190;
t143 = mrSges(1,1) + t184;
t133 = t143 * g(3);
t129 = qJ(2,1) * t231 + t292;
t128 = qJ(2,2) * t231 + t292;
t127 = qJ(2,3) * t231 + t292;
t117 = g(1) * t173 + g(2) * t176;
t116 = g(1) * t172 + g(2) * t175;
t115 = g(1) * t171 + g(2) * t174;
t72 = (t225 * t78 + t231 * t93) * t170;
t71 = (t223 * t77 + t231 * t92) * t169;
t70 = (t221 * t76 + t231 * t91) * t168;
t63 = t69 * t317;
t62 = t68 * t318;
t61 = t67 * t319;
t60 = (Ifges(3,3) * t332 + t75 * t350) * t153;
t59 = (Ifges(3,3) * t326 + t75 * t352) * t153;
t58 = (Ifges(3,3) * t334 + t74 * t354) * t152;
t57 = (Ifges(3,3) * t328 + t74 * t366) * t152;
t56 = (Ifges(3,3) * t336 + t73 * t367) * t151;
t55 = (Ifges(3,3) * t330 + t73 * t368) * t151;
t54 = (-0.2e1 * t320 + (t194 + t305) * t224 - 0.2e1 * t300 + t270) * t189 + (t196 * t224 + (t356 - t384) * t410 + (0.4e1 * t360 + t193 * t224 + 0.2e1 * (t206 * t224 - t402) * t218 + t388) * t201) * t202 + t320 + 0.2e1 * (-t303 - t381) * t224 + (t164 + mrSges(2,2)) * t312 + (t205 ^ 2) * m(3) + m(2) * t248 + qJ(2,1) * t405 + t254;
t53 = (-0.2e1 * t321 + (t194 + t306) * t222 - 0.2e1 * t301 + t270) * t189 + (t196 * t222 + (t356 - t385) * t410 + (0.4e1 * t361 + t193 * t222 + 0.2e1 * (t206 * t222 - t402) * t216 + t388) * t201) * t202 + t321 + 0.2e1 * (-t303 - t382) * t222 + (t163 + mrSges(2,2)) * t312 + (t204 ^ 2) * m(3) + m(2) * t247 + qJ(2,2) * t405 + t254;
t52 = (-0.2e1 * t322 + (t194 + t307) * t220 - 0.2e1 * t302 + t270) * t189 + (t196 * t220 + (t356 - t386) * t410 + (0.4e1 * t362 + t193 * t220 + 0.2e1 * (t206 * t220 - t402) * t214 + t388) * t201) * t202 + t322 + 0.2e1 * (-t303 - t383) * t220 + (t162 + mrSges(2,2)) * t312 + (t203 ^ 2) * m(3) + m(2) * t246 + qJ(2,3) * t405 + t254;
t39 = (t225 * t54 + t78 * t93) * t170;
t38 = (t223 * t53 + t77 * t92) * t169;
t37 = (t221 * t52 + t76 * t91) * t168;
t36 = (t51 * t344 + t78 * t351) * t170;
t35 = (t50 * t344 + t78 * t353) * t170;
t34 = (t49 * t346 + t77 * t355) * t169;
t33 = (t48 * t346 + t77 * t370) * t169;
t32 = (t47 * t348 + t76 * t372) * t168;
t31 = (t46 * t348 + t76 * t373) * t168;
t27 = t51 * t289 + (t75 * t332 + t54 * t350) * t153;
t26 = t50 * t289 + (t75 * t326 + t54 * t352) * t153;
t25 = t49 * t290 + (t74 * t334 + t53 * t354) * t152;
t24 = t48 * t290 + (t74 * t328 + t53 * t366) * t152;
t23 = t47 * t291 + (t73 * t336 + t52 * t367) * t151;
t22 = t46 * t291 + (t73 * t330 + t52 * t368) * t151;
t18 = t66 - t43 / 0.2e1 - t42 / 0.2e1 - t80 / 0.2e1;
t17 = t65 - t41 / 0.2e1 - t40 / 0.2e1 - t79 / 0.2e1;
t16 = t64 - t45 / 0.2e1 - t44 / 0.2e1 - t81 / 0.2e1;
t15 = (-pkin(3) * t369 + (-t64 + 0.2e1 * t30 + (-t160 - t395) * t69) * t69) * t170;
t14 = (-pkin(3) * t371 + (-t66 + 0.2e1 * t29 + (-t160 - t396) * t68) * t68) * t169;
t13 = (-pkin(3) * t374 + (-t65 + 0.2e1 * t28 + (-t160 - t397) * t67) * t67) * t168;
t12 = ((t16 * t412 + (qJ(2,1) * t404 - t249 * cos(t150) - 0.2e1 * t248 + t253) * t69 / 0.2e1 + (t114 + t411) * t30) * t69 + ((t90 * t188 * t156 - 0.2e1 * t16 * t159 + (-cos(qJ(3,1) + t241) * pkin(2) - t389) * t69) * t69 - (-t69 * t188 * sin(t150) / 0.2e1 + t90 * t114) * t153 * t90) * pkin(3)) * t170;
t11 = ((t18 * t412 + (qJ(2,2) * t404 - t249 * cos(t149) - 0.2e1 * t247 + t253) * t68 / 0.2e1 + (t113 + t411) * t29) * t68 + ((t89 * t187 * t155 - 0.2e1 * t18 * t158 + (-cos(t241 + qJ(3,2)) * pkin(2) - t390) * t68) * t68 - (-t68 * t187 * sin(t149) / 0.2e1 + t89 * t113) * t152 * t89) * pkin(3)) * t169;
t10 = ((t17 * t412 + (qJ(2,3) * t404 - t249 * cos(t148) - 0.2e1 * t246 + t253) * t67 / 0.2e1 + (t112 + t411) * t28) * t67 + ((t88 * t186 * t154 - 0.2e1 * t17 * t157 + (-cos(qJ(3,3) + t241) * pkin(2) - t391) * t67) * t67 - (-t67 * t186 * sin(t148) / 0.2e1 + t88 * t112) * t151 * t88) * pkin(3)) * t168;
t9 = -t75 * t15 + Ifges(3,3) * t286 + t69 * (-0.4e1 * (t360 + (t177 * t218 + t181) * t224 + t218 * t293 + t234) * t363 + (-0.2e1 * t271 + (t400 + (t305 + t402) * t357) * t224 - t280 + t19 * t164) * t202 + 0.2e1 * t296 + (t19 * t387 + t63) * t224 - t323 * t400 - Ifges(3,4) * t69) + (-mrSges(3,1) * t117 + t264 * mrSges(3,2)) * t159 + t156 * (t264 * mrSges(3,1) + mrSges(3,2) * t117);
t8 = -t74 * t14 + Ifges(3,3) * t287 + t68 * (-0.4e1 * (t361 + (t177 * t216 + t181) * t222 + t216 * t293 + t234) * t364 + (-0.2e1 * t272 + (t398 + (t306 + t402) * t358) * t222 - t281 + t21 * t163) * t202 + 0.2e1 * t297 + (t21 * t387 + t62) * t222 - t324 * t398 - Ifges(3,4) * t68) + (-mrSges(3,1) * t116 + t265 * mrSges(3,2)) * t158 + t155 * (t265 * mrSges(3,1) + mrSges(3,2) * t116);
t7 = -t73 * t13 + Ifges(3,3) * t288 + t67 * (-0.4e1 * (t362 + (t177 * t214 + t181) * t220 + t214 * t293 + t234) * t365 + (-0.2e1 * t273 + (t399 + (t307 + t402) * t359) * t220 - t282 + t20 * t162) * t202 + 0.2e1 * t298 + (t20 * t387 + t61) * t220 - t325 * t399 - Ifges(3,4) * t67) + (-mrSges(3,1) * t115 + t266 * mrSges(3,2)) * t157 + t154 * (t266 * mrSges(3,1) + mrSges(3,2) * t115);
t6 = -t78 * t15 - t69 * (t69 * t129 + 0.2e1 * (-t267 * t201 - t314 * t202) * t90) + (-g(3) * t219 + t120 * t225 - t12) * t231;
t5 = -t77 * t14 - t68 * (t68 * t128 + 0.2e1 * (-t268 * t201 - t315 * t202) * t89) + (-g(3) * t217 + t119 * t223 - t11) * t231;
t4 = -t76 * t13 - t67 * (t67 * t127 + 0.2e1 * (-t269 * t201 - t316 * t202) * t88) + (-g(3) * t215 + t118 * t221 - t10) * t231;
t3 = -t54 * t15 - t78 * t12 + t75 * t286 + 0.4e1 * (0.2e1 * t360 + (t182 + t317) * t224 + t218 * t294 - Ifges(3,4)) * t90 * t363 + (t271 + (((-qJ(2,1) / 0.4e1 + t242) * mrSges(3,1) + t233) * t90 + ((t183 + t381) * t407 + t295) * t69) * t224 + t280 / 0.2e1 - t218 * (-t139 * t90 + 0.2e1 * t308) / 0.4e1) * t90 * t406 - 0.4e1 * t90 * t296 - 0.2e1 * ((((-qJ(2,1) / 0.2e1 + t243) * mrSges(3,2) + t232) * t90 + t308) * t201 + t63) * t90 * t224 - t90 * (-t142 * t90 + t69 * t313) * t323 + 0.2e1 * t69 * (Ifges(3,4) * t90 + t129 * t30) + (-t257 * g(3) + (-t143 + t261) * t120) * t225 - t219 * (t261 * g(3) + t257 * t120 - t133);
t2 = -t53 * t14 - t77 * t11 + t74 * t287 + 0.4e1 * (0.2e1 * t361 + (t182 + t318) * t222 + t216 * t294 - Ifges(3,4)) * t89 * t364 + (t272 + (((-qJ(2,2) / 0.4e1 + t242) * mrSges(3,1) + t233) * t89 + ((t183 + t382) * t407 + t295) * t68) * t222 + t281 / 0.2e1 - t216 * (-t138 * t89 + 0.2e1 * t309) / 0.4e1) * t89 * t406 - 0.4e1 * t89 * t297 - 0.2e1 * ((((-qJ(2,2) / 0.2e1 + t243) * mrSges(3,2) + t232) * t89 + t309) * t201 + t62) * t89 * t222 - t89 * (-t141 * t89 + t68 * t313) * t324 + 0.2e1 * t68 * (Ifges(3,4) * t89 + t128 * t29) + (-t256 * g(3) + (-t143 + t262) * t119) * t223 - t217 * (t262 * g(3) + t256 * t119 - t133);
t1 = -t52 * t13 - t76 * t10 + t73 * t288 + 0.4e1 * (0.2e1 * t362 + (t182 + t319) * t220 + t214 * t294 - Ifges(3,4)) * t88 * t365 + (t273 + (((-qJ(2,3) / 0.4e1 + t242) * mrSges(3,1) + t233) * t88 + ((t183 + t383) * t407 + t295) * t67) * t220 + t282 / 0.2e1 - t214 * (-t137 * t88 + 0.2e1 * t310) / 0.4e1) * t88 * t406 - 0.4e1 * t88 * t298 - 0.2e1 * ((((-qJ(2,3) / 0.2e1 + t243) * mrSges(3,2) + t232) * t88 + t310) * t201 + t61) * t88 * t220 - t88 * (-t140 * t88 + t67 * t313) * t325 + 0.2e1 * t67 * (Ifges(3,4) * t88 + t127 * t28) + (-t255 * g(3) + (-t143 + t263) * t118) * t221 - t215 * (t263 * g(3) + t255 * t118 - t133);
t85 = [(-g(1) + t213) * m(4) + ((t27 * t351 + t36 * t375) * t213 + (t27 * t353 + t36 * t376) * t212 + (t225 * t27 + t36 * t93) * t211 + t3 * t351 + t6 * t375) * t170 + ((t176 * t212 * t60 + (t213 * t60 + t9) * t173) * t153 + (t175 * t212 * t58 + (t213 * t58 + t8) * t172) * t152 + (t174 * t212 * t56 + (t213 * t56 + t7) * t171) * t151) * t250 + ((t25 * t355 + t34 * t377) * t213 + (t25 * t370 + t34 * t378) * t212 + (t223 * t25 + t34 * t92) * t211 + t2 * t355 + t5 * t377) * t169 + ((t23 * t372 + t32 * t379) * t213 + (t23 * t373 + t32 * t380) * t212 + (t221 * t23 + t32 * t91) * t211 + t1 * t372 + t4 * t379) * t168; (-g(2) + t212) * m(4) + ((t26 * t351 + t35 * t375) * t213 + (t26 * t353 + t35 * t376) * t212 + (t225 * t26 + t35 * t93) * t211 + t3 * t353 + t6 * t376) * t170 + ((t173 * t213 * t59 + (t212 * t59 + t9) * t176) * t153 + (t172 * t213 * t57 + (t212 * t57 + t8) * t175) * t152 + (t171 * t213 * t55 + (t212 * t55 + t7) * t174) * t151) * t250 + ((t24 * t355 + t33 * t377) * t213 + (t24 * t370 + t33 * t378) * t212 + (t223 * t24 + t33 * t92) * t211 + t2 * t370 + t5 * t378) * t169 + ((t22 * t372 + t31 * t379) * t213 + (t22 * t373 + t31 * t380) * t212 + (t22 * t221 + t31 * t91) * t211 + t1 * t373 + t4 * t380) * t168; (-g(3) + t211) * m(4) + ((t225 * t39 + t72 * t93) * t211 + t225 * t3 + t93 * t6 + (t212 * t50 + t213 * t51) * t72 * t108 + ((t102 * t39 + t173 * t283) * t213 + (t101 * t39 + t176 * t283) * t212) * t153) * t170 + ((t223 * t38 + t71 * t92) * t211 + t223 * t2 + t92 * t5 + (t212 * t48 + t213 * t49) * t71 * t107 + ((t100 * t38 + t172 * t284) * t213 + (t175 * t284 + t38 * t99) * t212) * t152) * t169 + ((t221 * t37 + t70 * t91) * t211 + t221 * t1 + t91 * t4 + (t212 * t46 + t213 * t47) * t70 * t106 + ((t171 * t285 + t37 * t98) * t213 + (t174 * t285 + t37 * t97) * t212) * t151) * t168;];
tauX  = t85;
