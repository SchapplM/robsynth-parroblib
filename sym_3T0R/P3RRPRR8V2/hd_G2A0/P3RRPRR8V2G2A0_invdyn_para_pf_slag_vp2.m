% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR8V2G2A0
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2020-08-06 21:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:11:54
% EndTime: 2020-08-06 21:12:03
% DurationCPUTime: 8.19s
% Computational Cost: add. (31869->598), mult. (55041->930), div. (4491->12), fcn. (38313->50), ass. (0->370)
t214 = (-qJ(3,1) - pkin(5));
t143 = mrSges(3,1) * t214 + Ifges(3,5);
t209 = sin(pkin(7));
t210 = cos(pkin(7));
t401 = mrSges(3,2) * t214 + Ifges(3,6);
t418 = -t143 * t210 + t209 * t401;
t213 = (-qJ(3,2) - pkin(5));
t142 = mrSges(3,1) * t213 + Ifges(3,5);
t402 = mrSges(3,2) * t213 + Ifges(3,6);
t417 = -t142 * t210 + t209 * t402;
t212 = (-qJ(3,3) - pkin(5));
t141 = mrSges(3,1) * t212 + Ifges(3,5);
t403 = mrSges(3,2) * t212 + Ifges(3,6);
t416 = -t141 * t210 + t209 * t403;
t197 = t210 ^ 2;
t415 = 0.2e1 * t197;
t400 = 2 * pkin(1);
t320 = 2 * pkin(3);
t414 = pkin(2) / 0.2e1;
t413 = t400 / 0.2e1;
t203 = qJ(2,1) + pkin(7);
t168 = cos(t203);
t151 = pkin(3) * t168;
t237 = cos(qJ(2,1));
t189 = t237 * pkin(2);
t323 = t151 + t189;
t412 = pkin(1) + t323;
t202 = qJ(2,2) + pkin(7);
t167 = cos(t202);
t150 = pkin(3) * t167;
t235 = cos(qJ(2,2));
t188 = t235 * pkin(2);
t324 = t150 + t188;
t411 = pkin(1) + t324;
t201 = qJ(2,3) + pkin(7);
t166 = cos(t201);
t149 = pkin(3) * t166;
t233 = cos(qJ(2,3));
t187 = t233 * pkin(2);
t325 = t149 + t187;
t410 = pkin(1) + t325;
t169 = t209 * mrSges(3,1);
t218 = Ifges(3,2) - Ifges(3,1);
t101 = (-pkin(2) * mrSges(3,2) - t209 * t218) * t210 - pkin(2) * t169 + Ifges(3,4) * t415 + Ifges(2,4) - Ifges(3,4);
t171 = (t212 * m(3));
t220 = (mrSges(2,3) + mrSges(3,3));
t292 = m(2) * pkin(5) - mrSges(1,2) + t220;
t406 = -t292 + t171;
t172 = (t213 * m(3));
t405 = -t292 + t172;
t173 = (t214 * m(3));
t404 = -t292 + t173;
t122 = 0.1e1 / t325;
t221 = legFrame(3,2);
t177 = sin(t221);
t180 = cos(t221);
t242 = xDP(2);
t243 = xDP(1);
t85 = (t177 * t243 + t180 * t242) * t122;
t82 = t85 ^ 2;
t123 = 0.1e1 / t324;
t222 = legFrame(2,2);
t178 = sin(t222);
t181 = cos(t222);
t86 = (t178 * t243 + t181 * t242) * t123;
t83 = t86 ^ 2;
t124 = 0.1e1 / t323;
t223 = legFrame(1,2);
t179 = sin(t223);
t182 = cos(t223);
t87 = (t179 * t243 + t182 * t242) * t124;
t84 = t87 ^ 2;
t399 = 0.2e1 * mrSges(3,1);
t398 = 2 * mrSges(3,3);
t262 = pkin(3) ^ 2;
t370 = pkin(3) * t210;
t161 = pkin(2) * t370;
t263 = pkin(2) ^ 2;
t211 = t263 / 0.2e1;
t322 = t161 + t211;
t397 = -0.4e1 * pkin(1) * (t262 / 0.2e1 + t322);
t396 = 0.2e1 * t161;
t205 = t233 ^ 2;
t395 = 0.2e1 * t205;
t206 = t235 ^ 2;
t394 = 0.2e1 * t206;
t207 = t237 ^ 2;
t393 = 0.2e1 * t207;
t392 = 0.2e1 * t210;
t227 = sin(qJ(2,3));
t391 = -0.2e1 * t227;
t229 = sin(qJ(2,2));
t390 = -0.2e1 * t229;
t231 = sin(qJ(2,1));
t389 = -0.2e1 * t231;
t388 = -0.4e1 * t233;
t387 = -0.4e1 * t235;
t386 = -0.4e1 * t237;
t385 = -4 * pkin(5) - 4 * pkin(6);
t249 = m(3) * pkin(2);
t384 = mrSges(2,2) * pkin(1);
t154 = pkin(2) + t370;
t331 = t209 * t227;
t279 = pkin(3) * t331 - t154 * t233;
t108 = 0.1e1 / t279;
t194 = pkin(6) - t212;
t174 = 1 / t194;
t234 = cos(qJ(1,3));
t241 = xDP(3);
t228 = sin(qJ(1,3));
t371 = pkin(3) * t209;
t308 = t228 * t371;
t311 = t180 * t371;
t337 = t177 * t228;
t340 = t154 * t180;
t76 = (-t154 * t337 + t311) * t233 + (t177 * t308 + t340) * t227;
t314 = t177 * t371;
t334 = t180 * t228;
t343 = t154 * t177;
t79 = (t154 * t334 + t314) * t233 + t227 * (-t180 * t308 + t343);
t64 = (t234 * t241 - (t242 * t76 + t243 * t79) * t108) * t174;
t62 = pkin(1) * t64;
t330 = t209 * t229;
t278 = pkin(3) * t330 - t154 * t235;
t109 = 0.1e1 / t278;
t195 = pkin(6) - t213;
t175 = 1 / t195;
t236 = cos(qJ(1,2));
t230 = sin(qJ(1,2));
t307 = t230 * t371;
t310 = t181 * t371;
t336 = t178 * t230;
t339 = t154 * t181;
t77 = (-t154 * t336 + t310) * t235 + (t178 * t307 + t339) * t229;
t313 = t178 * t371;
t333 = t181 * t230;
t342 = t154 * t178;
t80 = (t154 * t333 + t313) * t235 + t229 * (-t181 * t307 + t342);
t65 = (t236 * t241 - (t242 * t77 + t243 * t80) * t109) * t175;
t63 = pkin(1) * t65;
t290 = pkin(1) * t228 - t194 * t234;
t305 = pkin(3) * (t197 - 0.1e1);
t383 = pkin(3) * (t228 * t305 + t290 * t331);
t289 = pkin(1) * t230 - t195 * t236;
t382 = pkin(3) * (t230 * t305 + t289 * t330);
t232 = sin(qJ(1,1));
t196 = pkin(6) - t214;
t238 = cos(qJ(1,1));
t288 = pkin(1) * t232 - t196 * t238;
t329 = t209 * t231;
t381 = pkin(3) * (t232 * t305 + t288 * t329);
t380 = pkin(5) * mrSges(2,2);
t277 = pkin(3) * t329 - t154 * t237;
t110 = 0.1e1 / t277;
t176 = 1 / t196;
t306 = t232 * t371;
t309 = t182 * t371;
t335 = t179 * t232;
t338 = t154 * t182;
t78 = (-t154 * t335 + t309) * t237 + (t179 * t306 + t338) * t231;
t312 = t179 * t371;
t332 = t182 * t232;
t341 = t154 * t179;
t81 = (t154 * t332 + t312) * t237 + t231 * (-t182 * t306 + t341);
t66 = (t238 * t241 - (t242 * t78 + t243 * t81) * t110) * t176;
t61 = t66 * pkin(1);
t379 = t122 / 0.2e1;
t378 = t123 / 0.2e1;
t377 = t124 / 0.2e1;
t376 = m(3) * t263;
t183 = mrSges(2,1) + t249;
t375 = pkin(1) * t183;
t374 = pkin(1) * t227;
t373 = pkin(1) * t229;
t372 = pkin(1) * t231;
t369 = pkin(3) * t263;
t368 = t262 * pkin(2);
t367 = t64 * t85;
t366 = t65 * t86;
t365 = t66 * t87;
t170 = mrSges(3,1) * t210;
t364 = mrSges(3,2) * t209;
t363 = mrSges(3,2) * t210;
t361 = Ifges(3,4) * t209;
t328 = t218 * t197;
t91 = 0.2e1 * t328 + (pkin(2) * t399 + 0.4e1 * t361) * t210 + t376 - 0.2e1 * pkin(2) * t364 + Ifges(2,2) - Ifges(2,1) - t218;
t360 = t227 * t91;
t359 = t229 * t91;
t358 = t231 * t91;
t357 = pkin(5) * mrSges(2,1) - Ifges(2,5);
t356 = -t183 * pkin(5) + Ifges(2,5);
t148 = t169 + mrSges(2,2);
t355 = t101 * t205;
t354 = t101 * t206;
t353 = t101 * t207;
t352 = t108 * t174;
t351 = t109 * t175;
t350 = t110 * t176;
t349 = t122 * t177;
t348 = t122 * t180;
t347 = t123 * t178;
t346 = t123 * t181;
t345 = t124 * t179;
t344 = t124 * t182;
t327 = -0.2e1 * t249;
t326 = pkin(2) * t320;
t321 = t262 + t263;
t319 = -0.2e1 * t369;
t318 = -0.2e1 * t148 * pkin(1);
t317 = -0.2e1 * t368;
t105 = t308 * t391 + t290;
t112 = (t197 - 0.1e1 / 0.2e1) * t262 + t322;
t132 = -t371 + t374;
t162 = pkin(1) * t371;
t271 = t262 * t415 - t262 + t263 + t396;
t98 = t227 * t271 + t162;
t55 = (-t112 * t337 + t154 * t311) * t395 + (-t105 * t343 + t180 * t98) * t233 + t177 * t383 + t132 * t340;
t40 = t55 * t242 * t352;
t56 = (t112 * t334 + t154 * t314) * t395 + (t105 * t340 + t177 * t98) * t233 - t180 * t383 + t132 * t343;
t41 = t56 * t243 * t352;
t102 = t228 * t194 + t234 * t410;
t92 = t102 * t174 * t241;
t28 = -t41 - t40 + t92;
t106 = t307 * t390 + t289;
t133 = -t371 + t373;
t99 = t229 * t271 + t162;
t57 = (-t112 * t336 + t154 * t310) * t394 + (-t106 * t342 + t181 * t99) * t235 + t178 * t382 + t133 * t339;
t42 = t57 * t242 * t351;
t58 = (t112 * t333 + t154 * t313) * t394 + (t106 * t339 + t178 * t99) * t235 - t181 * t382 + t133 * t342;
t43 = t58 * t243 * t351;
t103 = t230 * t195 + t236 * t411;
t93 = t103 * t175 * t241;
t29 = -t43 - t42 + t93;
t100 = t231 * t271 + t162;
t107 = t306 * t389 + t288;
t134 = -t371 + t372;
t59 = (-t112 * t335 + t154 * t309) * t393 + (t100 * t182 - t107 * t341) * t237 + t179 * t381 + t134 * t338;
t44 = t59 * t242 * t350;
t60 = (t112 * t332 + t154 * t312) * t393 + (t100 * t179 + t107 * t338) * t237 - t182 * t381 + t134 * t341;
t45 = t60 * t243 * t350;
t104 = t232 * t196 + t238 * t412;
t94 = t104 * t176 * t241;
t30 = -t45 - t44 + t94;
t295 = -m(3) * qJ(3,3) - mrSges(3,3);
t298 = Ifges(2,6) - t380;
t73 = (pkin(2) * t295 + t356 - t416) * t227 + t233 * (t141 * t209 + t210 * t403 + t298);
t304 = t73 * t352;
t296 = -m(3) * qJ(3,2) - mrSges(3,3);
t74 = (pkin(2) * t296 + t356 - t417) * t229 + t235 * (t142 * t209 + t210 * t402 + t298);
t303 = t74 * t351;
t297 = -m(3) * qJ(3,1) - mrSges(3,3);
t75 = (pkin(2) * t297 + t356 - t418) * t231 + t237 * (t143 * t209 + t210 * t401 + t298);
t302 = t75 * t350;
t126 = pkin(2) * t227 + pkin(3) * sin(t201);
t301 = t122 * t126 * t82;
t127 = pkin(2) * t229 + pkin(3) * sin(t202);
t300 = t123 * t127 * t83;
t128 = pkin(2) * t231 + pkin(3) * sin(t203);
t299 = t124 * t128 * t84;
t294 = -(2 * pkin(3) * t262) - 0.4e1 * t369;
t293 = t249 - t364;
t291 = pkin(1) * t399 * t210 + (mrSges(2,1) + t293) * t400;
t114 = -t170 - t293;
t287 = -t364 + t170;
t286 = t169 + t363;
t118 = g(1) * t180 - g(2) * t177;
t285 = -g(3) * t234 - t118 * t228;
t119 = g(1) * t181 - g(2) * t178;
t284 = -g(3) * t236 - t119 * t230;
t120 = g(1) * t182 - g(2) * t179;
t283 = -g(3) * t238 - t120 * t232;
t113 = -mrSges(2,1) + t114;
t129 = t148 + t363;
t282 = -t113 * t233 - t129 * t227;
t281 = -t113 * t235 - t129 * t229;
t280 = -t113 * t237 - t129 * t231;
t260 = pkin(5) ^ 2;
t264 = pkin(1) ^ 2;
t276 = -(2 * t260) - (2 * t264) - t321 + ((-4 * pkin(5) - 2 * pkin(6)) * pkin(6));
t275 = 0.2e1 * t101;
t274 = t286 * t233;
t273 = t286 * t235;
t272 = t286 * t237;
t269 = (mrSges(2,2) + t286) * t400;
t267 = (m(2) * t260) + (2 * t220 * pkin(5)) + (m(2) + m(3)) * t264 + Ifges(2,1) + Ifges(3,2) + Ifges(1,3) - t328;
t259 = 0.2e1 * qJ(2,1);
t258 = 0.2e1 * qJ(2,2);
t257 = 0.2e1 * qJ(2,3);
t256 = -pkin(5) / 0.2e1;
t255 = 0.2e1 * pkin(7);
t250 = m(3) * pkin(1);
t248 = m(3) * pkin(5);
t247 = Ifges(3,5) / 0.2e1;
t246 = Ifges(3,6) / 0.2e1;
t226 = xDDP(1);
t225 = xDDP(2);
t224 = xDDP(3);
t200 = t257 + pkin(7);
t199 = pkin(7) + t259;
t198 = pkin(7) + t258;
t186 = -qJ(3,1) / 0.2e1 + t256;
t185 = -qJ(3,2) / 0.2e1 + t256;
t184 = -qJ(3,3) / 0.2e1 + t256;
t160 = pkin(2) * t211 + t368;
t158 = -t380 / 0.2e1 + Ifges(2,6) / 0.2e1;
t157 = 0.2e1 * t203;
t156 = 0.2e1 * t202;
t155 = 0.2e1 * t201;
t147 = (m(2) * pkin(1)) + mrSges(1,1) + t250;
t146 = t248 - t297;
t145 = t248 - t296;
t144 = t248 - t295;
t137 = t147 * g(3);
t135 = t396 + t321;
t117 = g(1) * t179 + g(2) * t182;
t116 = g(1) * t178 + g(2) * t181;
t115 = g(1) * t177 + g(2) * t180;
t111 = 0.2e1 * pkin(2) * t287 + Ifges(2,3) + Ifges(3,3) + t376;
t90 = t114 * t237 + t231 * t286 - t250;
t89 = t114 * t235 + t229 * t286 - t250;
t88 = t114 * t233 + t227 * t286 - t250;
t72 = (m(3) * t104 + t238 * t90) * t176;
t71 = (m(3) * t103 + t236 * t89) * t175;
t70 = (m(3) * t102 + t234 * t88) * t174;
t69 = t91 * t207 + (t231 * t275 + t291) * t237 + (-mrSges(3,2) * t372 - t361) * t392 + t231 * t318 + (t214 ^ 2) * m(3) + qJ(3,1) * t398 + t267;
t68 = t91 * t206 + (t229 * t275 + t291) * t235 + (-mrSges(3,2) * t373 - t361) * t392 + t229 * t318 + (t213 ^ 2) * m(3) + qJ(3,2) * t398 + t267;
t67 = t91 * t205 + (t227 * t275 + t291) * t233 + (-mrSges(3,2) * t374 - t361) * t392 + t227 * t318 + (t212 ^ 2) * m(3) + qJ(3,3) * t398 + t267;
t54 = t66 * t375;
t53 = t65 * t375;
t52 = t64 * t375;
t51 = t111 * t345 - t302 * t81;
t50 = t111 * t347 - t303 * t80;
t49 = t111 * t349 - t304 * t79;
t48 = t111 * t344 - t302 * t78;
t47 = t111 * t346 - t303 * t77;
t46 = t111 * t348 - t304 * t76;
t39 = (t104 * t90 + t238 * t69) * t176;
t38 = (t103 * t89 + t236 * t68) * t175;
t37 = (t102 * t88 + t234 * t67) * t174;
t36 = (m(3) * t60 + t81 * t90) * t350;
t35 = (m(3) * t58 + t80 * t89) * t351;
t34 = (m(3) * t56 + t79 * t88) * t352;
t33 = (m(3) * t59 + t78 * t90) * t350;
t32 = (m(3) * t57 + t77 * t89) * t351;
t31 = (m(3) * t55 + t76 * t88) * t352;
t27 = t75 * t345 - (t60 * t90 + t69 * t81) * t350;
t26 = t74 * t347 - (t58 * t89 + t68 * t80) * t351;
t25 = t73 * t349 - (t56 * t88 + t67 * t79) * t352;
t24 = t75 * t344 - (t59 * t90 + t69 * t78) * t350;
t23 = t74 * t346 - (t57 * t89 + t68 * t77) * t351;
t22 = t73 * t348 - (t55 * t88 + t67 * t76) * t352;
t18 = t63 + t43 / 0.2e1 + t42 / 0.2e1 - t93 / 0.2e1;
t17 = t62 + t41 / 0.2e1 + t40 / 0.2e1 - t92 / 0.2e1;
t16 = t61 + t45 / 0.2e1 + t44 / 0.2e1 - t94 / 0.2e1;
t15 = (-t84 * t135 / (t189 + (t210 * t237 - t329) * pkin(3)) + (t277 * t66 + 0.2e1 * t30 - t61) * t66) * t176;
t14 = (-t83 * t135 / (t188 + (t210 * t235 - t330) * pkin(3)) + (t278 * t65 + 0.2e1 * t29 - t63) * t65) * t175;
t13 = (-t82 * t135 / (t187 + (t210 * t233 - t331) * pkin(3)) + (t279 * t64 + 0.2e1 * t28 - t62) * t64) * t174;
t12 = ((cos(qJ(2,1) - pkin(7)) * t319 + cos(qJ(2,1) + t255) * t317 + t294 * t168 + t160 * t386 + t397) * t84 * t377 + (-0.2e1 * t16 * t151 + (-t262 * cos(t157) - t263 * cos(t259) + t276 + (t385 - 0.2e1 * qJ(3,1)) * qJ(3,1)) * t66 / 0.2e1 + (t16 * t386 + (-cos(t199) - t210) * t66 * t320) * t414 + (t128 + (sin(t199) * t326 + sin(t157) * t262 + sin(t259) * t263) * t377) * t87 * t196 + (t412 + t413) * t30) * t66) * t176;
t11 = ((cos(qJ(2,2) - pkin(7)) * t319 + cos(t255 + qJ(2,2)) * t317 + t294 * t167 + t160 * t387 + t397) * t83 * t378 + (-0.2e1 * t18 * t150 + (-cos(t156) * t262 - cos(t258) * t263 + t276 + (t385 - 0.2e1 * qJ(3,2)) * qJ(3,2)) * t65 / 0.2e1 + (t18 * t387 + (-cos(t198) - t210) * t65 * t320) * t414 + (t127 + (sin(t198) * t326 + sin(t156) * t262 + sin(t258) * t263) * t378) * t86 * t195 + (t411 + t413) * t29) * t65) * t175;
t10 = ((cos(qJ(2,3) - pkin(7)) * t319 + cos(t255 + qJ(2,3)) * t317 + t294 * t166 + t160 * t388 + t397) * t82 * t379 + (-0.2e1 * t17 * t149 + (-cos(t155) * t262 - cos(t257) * t263 + t276 + (t385 - 0.2e1 * qJ(3,3)) * qJ(3,3)) * t64 / 0.2e1 + (t17 * t388 + (-cos(t200) - t210) * t64 * t320) * t414 + (t126 + (sin(t200) * t326 + sin(t155) * t262 + sin(t257) * t263) * t379) * t85 * t194 + (t410 + t413) * t28) * t64) * t174;
t9 = -t75 * t15 + t111 * t299 + ((t30 * t327 + t54) * t231 + (t231 * t287 + t272) * (t61 + 0.2e1 * t45 + 0.2e1 * t44 - 0.2e1 * t94) + (-0.2e1 * t353 + (t358 + t384) * t237 + t101) * t66) * t66 + (t283 * t113 + t117 * t129) * t231 - t237 * (-t113 * t117 + t283 * t129);
t8 = -t74 * t14 + t111 * t300 + ((t29 * t327 + t53) * t229 + (t229 * t287 + t273) * (t63 + 0.2e1 * t43 + 0.2e1 * t42 - 0.2e1 * t93) + (-0.2e1 * t354 + (t359 + t384) * t235 + t101) * t65) * t65 + (t284 * t113 + t116 * t129) * t229 - t235 * (-t113 * t116 + t284 * t129);
t7 = -t73 * t13 + t111 * t301 + ((t28 * t327 + t52) * t227 + (t227 * t287 + t274) * (t62 + 0.2e1 * t41 + 0.2e1 * t40 - 0.2e1 * t92) + (-0.2e1 * t355 + (t360 + t384) * t233 + t101) * t64) * t64 + (t285 * t113 + t115 * t129) * t227 - t233 * (-t113 * t115 + t285 * t129);
t6 = -t90 * t15 - t66 ^ 2 * (-t173 + mrSges(3,3)) + (-g(3) * t232 + t120 * t238 - t12) * m(3) + 0.2e1 * (-t114 * t231 + t272) * t365;
t5 = -t89 * t14 - t65 ^ 2 * (-t172 + mrSges(3,3)) + (-g(3) * t230 + t119 * t236 - t11) * m(3) + 0.2e1 * (-t114 * t229 + t273) * t366;
t4 = -t88 * t13 - t64 ^ 2 * (-t171 + mrSges(3,3)) + (-g(3) * t228 + t118 * t234 - t10) * m(3) + 0.2e1 * (-t114 * t227 + t274) * t367;
t3 = -t69 * t15 + t75 * t299 - t90 * t12 + 0.4e1 * t353 * t365 - ((pkin(2) * t146 + t357 + t418) * t87 + (t269 + 0.2e1 * t358) * t66) * t87 * t237 + (((mrSges(3,2) * t186 + t246) * t87 + mrSges(3,1) * t61) * t210 + ((mrSges(3,1) * t186 + t247) * t87 - mrSges(3,2) * t61) * t209 + t158 * t87 + t54) * t87 * t389 + 0.2e1 * t66 * (-t101 * t87 + t146 * t30) + (t404 * g(3) + (-t147 - t280) * t120) * t238 + (t280 * g(3) + t120 * t404 + t137) * t232;
t2 = -t68 * t14 + t74 * t300 - t89 * t11 + 0.4e1 * t354 * t366 - ((pkin(2) * t145 + t357 + t417) * t86 + (t269 + 0.2e1 * t359) * t65) * t86 * t235 + (((mrSges(3,2) * t185 + t246) * t86 + mrSges(3,1) * t63) * t210 + ((mrSges(3,1) * t185 + t247) * t86 - mrSges(3,2) * t63) * t209 + t158 * t86 + t53) * t86 * t390 + 0.2e1 * t65 * (-t101 * t86 + t145 * t29) + (t405 * g(3) + (-t147 - t281) * t119) * t236 + (t281 * g(3) + t119 * t405 + t137) * t230;
t1 = -t67 * t13 + t73 * t301 - t88 * t10 + 0.4e1 * t355 * t367 - ((pkin(2) * t144 + t357 + t416) * t85 + (t269 + 0.2e1 * t360) * t64) * t85 * t233 + (((mrSges(3,2) * t184 + t246) * t85 + mrSges(3,1) * t62) * t210 + ((mrSges(3,1) * t184 + t247) * t85 - mrSges(3,2) * t62) * t209 + t158 * t85 + t52) * t85 * t391 + 0.2e1 * t64 * (-t101 * t85 + t144 * t28) + (t406 * g(3) + (-t147 - t282) * t118) * t234 + (t282 * g(3) + t118 * t406 + t137) * t228;
t19 = [t7 * t349 + t8 * t347 + t9 * t345 - g(1) * m(4) + (t344 * t51 + t346 * t50 + t348 * t49) * t225 + (t345 * t51 + t347 * t50 + t349 * t49 + m(4)) * t226 + ((-t104 * t36 + t238 * t27) * t224 - ((t27 * t81 - t36 * t60) * t226 + (t27 * t78 - t36 * t59) * t225 + t81 * t3 + t60 * t6) * t110) * t176 + ((-t103 * t35 + t236 * t26) * t224 - ((t26 * t80 - t35 * t58) * t226 + (t26 * t77 - t35 * t57) * t225 + t80 * t2 + t58 * t5) * t109) * t175 + ((-t102 * t34 + t234 * t25) * t224 - ((t25 * t79 - t34 * t56) * t226 + (t25 * t76 - t34 * t55) * t225 + t79 * t1 + t56 * t4) * t108) * t174; t7 * t348 + t8 * t346 + t9 * t344 - g(2) * m(4) + (t345 * t48 + t347 * t47 + t349 * t46) * t226 + (t344 * t48 + t346 * t47 + t348 * t46 + m(4)) * t225 + ((-t104 * t33 + t238 * t24) * t224 - ((t24 * t81 - t33 * t60) * t226 + (t24 * t78 - t33 * t59) * t225 + t78 * t3 + t59 * t6) * t110) * t176 + ((-t103 * t32 + t23 * t236) * t224 - ((t23 * t80 - t32 * t58) * t226 + (t23 * t77 - t32 * t57) * t225 + t77 * t2 + t57 * t5) * t109) * t175 + ((-t102 * t31 + t22 * t234) * t224 - ((t22 * t79 - t31 * t56) * t226 + (t22 * t76 - t31 * t55) * t225 + t76 * t1 + t55 * t4) * t108) * t174; (-g(3) + t224) * m(4) + ((t72 * t224 + t6) * t104 + (t39 * t224 + t3 + (t179 * t226 + t182 * t225) * t75 * t124) * t238 - ((t39 * t81 + t60 * t72) * t226 + (t39 * t78 + t59 * t72) * t225) * t110) * t176 + ((t71 * t224 + t5) * t103 + (t38 * t224 + t2 + (t178 * t226 + t181 * t225) * t74 * t123) * t236 - ((t38 * t80 + t58 * t71) * t226 + (t38 * t77 + t57 * t71) * t225) * t109) * t175 + ((t70 * t224 + t4) * t102 + (t37 * t224 + t1 + (t177 * t226 + t180 * t225) * t73 * t122) * t234 - ((t37 * t79 + t56 * t70) * t226 + (t37 * t76 + t55 * t70) * t225) * t108) * t174;];
tauX  = t19;
