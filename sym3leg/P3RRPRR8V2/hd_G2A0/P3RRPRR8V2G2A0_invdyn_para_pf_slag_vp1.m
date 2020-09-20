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
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
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

function tauX = P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:11:35
% EndTime: 2020-08-06 21:11:43
% DurationCPUTime: 7.75s
% Computational Cost: add. (32277->624), mult. (55506->965), div. (4569->12), fcn. (37437->50), ass. (0->401)
t453 = pkin(1) / 0.2e1;
t452 = pkin(2) / 0.2e1;
t263 = 2 * qJ(2,1);
t210 = t263 + pkin(7);
t174 = sin(t210);
t180 = cos(t210);
t451 = rSges(3,1) * t174 + rSges(3,2) * t180;
t261 = 2 * qJ(2,3);
t207 = t261 + pkin(7);
t171 = sin(t207);
t177 = cos(t207);
t450 = rSges(3,1) * t171 + rSges(3,2) * t177;
t262 = 2 * qJ(2,2);
t206 = pkin(7) + t262;
t170 = sin(t206);
t176 = cos(t206);
t449 = rSges(3,1) * t170 + rSges(3,2) * t176;
t448 = -0.2e1 * pkin(1);
t447 = 0.2e1 * pkin(1);
t446 = 0.2e1 * pkin(2);
t208 = qJ(2,3) + pkin(7);
t178 = cos(t208);
t154 = pkin(3) * t178;
t237 = cos(qJ(2,3));
t192 = t237 * pkin(2);
t348 = t154 + t192;
t121 = 0.1e1 / t348;
t225 = legFrame(3,2);
t185 = sin(t225);
t188 = cos(t225);
t251 = xDP(2);
t252 = xDP(1);
t82 = (t185 * t252 + t188 * t251) * t121;
t79 = t82 ^ 2;
t209 = qJ(2,2) + pkin(7);
t179 = cos(t209);
t155 = pkin(3) * t179;
t239 = cos(qJ(2,2));
t193 = t239 * pkin(2);
t347 = t155 + t193;
t122 = 0.1e1 / t347;
t226 = legFrame(2,2);
t186 = sin(t226);
t189 = cos(t226);
t83 = (t186 * t252 + t189 * t251) * t122;
t80 = t83 ^ 2;
t211 = qJ(2,1) + pkin(7);
t181 = cos(t211);
t156 = pkin(3) * t181;
t241 = cos(qJ(2,1));
t194 = t241 * pkin(2);
t346 = t156 + t194;
t123 = 0.1e1 / t346;
t227 = legFrame(1,2);
t187 = sin(t227);
t190 = cos(t227);
t84 = (t187 * t252 + t190 * t251) * t123;
t81 = t84 ^ 2;
t223 = cos(pkin(7));
t408 = pkin(3) * t223;
t157 = pkin(2) + t408;
t222 = sin(pkin(7));
t231 = sin(qJ(2,3));
t353 = t222 * t231;
t282 = pkin(3) * t353 - t157 * t237;
t107 = 0.1e1 / t282;
t397 = pkin(5) + qJ(3,3);
t202 = pkin(6) + t397;
t182 = 0.1e1 / t202;
t386 = t107 * t182;
t232 = sin(qJ(1,3));
t238 = cos(qJ(1,3));
t299 = pkin(1) * t232 - t202 * t238;
t409 = pkin(3) * t222;
t314 = t232 * t409;
t437 = -0.2e1 * t231;
t101 = t314 * t437 + t299;
t205 = t223 ^ 2;
t270 = pkin(3) ^ 2;
t167 = pkin(2) * t408;
t271 = pkin(2) ^ 2;
t224 = t271 / 0.2e1;
t345 = t167 + t224;
t111 = (t205 - 0.1e1 / 0.2e1) * t270 + t345;
t128 = pkin(1) * t231 - t409;
t317 = t188 * t409;
t359 = t185 * t232;
t368 = t157 * t188;
t371 = t157 * t185;
t311 = pkin(3) * (t205 - 0.1e1);
t426 = pkin(3) * (t232 * t311 + t299 * t353);
t440 = 0.2e1 * t237 ^ 2;
t168 = pkin(1) * t409;
t441 = 0.2e1 * t167;
t277 = 0.2e1 * t205 * t270 - t270 + t271 + t441;
t95 = t277 * t231 + t168;
t55 = (-t111 * t359 + t157 * t317) * t440 + (-t101 * t371 + t188 * t95) * t237 + t185 * t426 + t128 * t368;
t46 = t55 * t251 * t386;
t320 = t185 * t409;
t356 = t188 * t232;
t56 = (t111 * t356 + t157 * t320) * t440 + (t101 * t368 + t185 * t95) * t237 - t188 * t426 + t128 * t371;
t47 = t56 * t252 * t386;
t250 = xDP(3);
t396 = t192 + pkin(1);
t98 = t232 * t202 + (t396 + t154) * t238;
t88 = t98 * t182 * t250;
t28 = -t47 - t46 + t88;
t445 = 0.2e1 * t28;
t233 = sin(qJ(2,2));
t352 = t222 * t233;
t281 = pkin(3) * t352 - t157 * t239;
t108 = 0.1e1 / t281;
t398 = pkin(5) + qJ(3,2);
t203 = pkin(6) + t398;
t183 = 0.1e1 / t203;
t385 = t108 * t183;
t234 = sin(qJ(1,2));
t240 = cos(qJ(1,2));
t298 = pkin(1) * t234 - t203 * t240;
t313 = t234 * t409;
t436 = -0.2e1 * t233;
t102 = t313 * t436 + t298;
t129 = pkin(1) * t233 - t409;
t316 = t189 * t409;
t358 = t186 * t234;
t367 = t157 * t189;
t370 = t157 * t186;
t425 = pkin(3) * (t234 * t311 + t298 * t352);
t439 = 0.2e1 * t239 ^ 2;
t96 = t277 * t233 + t168;
t57 = (-t111 * t358 + t157 * t316) * t439 + (-t102 * t370 + t189 * t96) * t239 + t186 * t425 + t129 * t367;
t48 = t57 * t251 * t385;
t319 = t186 * t409;
t355 = t189 * t234;
t58 = (t111 * t355 + t157 * t319) * t439 + (t102 * t367 + t186 * t96) * t239 - t189 * t425 + t129 * t370;
t49 = t58 * t252 * t385;
t395 = t193 + pkin(1);
t99 = t234 * t203 + (t395 + t155) * t240;
t89 = t99 * t183 * t250;
t29 = -t49 - t48 + t89;
t444 = 0.2e1 * t29;
t235 = sin(qJ(2,1));
t351 = t222 * t235;
t280 = pkin(3) * t351 - t157 * t241;
t109 = 0.1e1 / t280;
t399 = pkin(5) + qJ(3,1);
t204 = pkin(6) + t399;
t184 = 0.1e1 / t204;
t384 = t109 * t184;
t236 = sin(qJ(1,1));
t242 = cos(qJ(1,1));
t297 = pkin(1) * t236 - t204 * t242;
t312 = t236 * t409;
t435 = -0.2e1 * t235;
t103 = t312 * t435 + t297;
t130 = pkin(1) * t235 - t409;
t315 = t190 * t409;
t357 = t187 * t236;
t366 = t157 * t190;
t369 = t157 * t187;
t424 = pkin(3) * (t236 * t311 + t297 * t351);
t438 = 0.2e1 * t241 ^ 2;
t97 = t277 * t235 + t168;
t59 = (-t111 * t357 + t157 * t315) * t438 + (-t103 * t369 + t190 * t97) * t241 + t187 * t424 + t130 * t366;
t50 = t59 * t251 * t384;
t318 = t187 * t409;
t354 = t190 * t236;
t60 = (t111 * t354 + t157 * t318) * t438 + (t103 * t366 + t187 * t97) * t241 - t190 * t424 + t130 * t369;
t51 = t60 * t252 * t384;
t394 = t194 + pkin(1);
t100 = t236 * t204 + (t394 + t156) * t242;
t90 = t100 * t184 * t250;
t30 = -t51 - t50 + t90;
t443 = 0.2e1 * t30;
t442 = -0.4e1 * pkin(1) * (t270 / 0.2e1 + t345);
t434 = -0.4e1 * t237;
t433 = -0.4e1 * t239;
t432 = -0.4e1 * t241;
t431 = -4 * pkin(5) - 4 * pkin(6);
t430 = m(3) * pkin(1);
t254 = pkin(2) * m(3);
t248 = m(2) * rSges(2,1);
t429 = (m(2) * rSges(2,2));
t428 = m(3) * rSges(3,2);
t427 = rSges(2,2) * pkin(1);
t423 = t121 / 0.2e1;
t422 = t122 / 0.2e1;
t421 = t123 / 0.2e1;
t265 = rSges(2,2) ^ 2;
t267 = rSges(2,1) ^ 2;
t124 = m(3) * t271 + (-t265 + t267) * m(2) + Icges(2,2) - Icges(2,1);
t420 = t124 / 0.2e1;
t264 = rSges(3,2) ^ 2;
t266 = rSges(3,1) ^ 2;
t134 = m(3) * (-t264 + t266) - Icges(3,1) + Icges(3,2);
t419 = t134 / 0.2e1;
t243 = rSges(2,3) + pkin(5);
t417 = m(2) * t243;
t172 = sin(t208);
t104 = -rSges(3,1) * t178 + rSges(3,2) * t172 - t396;
t416 = m(3) * t104;
t173 = sin(t209);
t105 = -rSges(3,1) * t179 + rSges(3,2) * t173 - t395;
t415 = m(3) * t105;
t175 = sin(t211);
t106 = -rSges(3,1) * t181 + rSges(3,2) * t175 - t394;
t414 = m(3) * t106;
t161 = t248 + t254;
t413 = pkin(1) * t161;
t412 = pkin(2) * t233;
t411 = pkin(2) * t235;
t410 = pkin(2) * t270;
t407 = pkin(3) * t271;
t199 = rSges(3,3) + t397;
t406 = t199 * m(3);
t200 = rSges(3,3) + t398;
t405 = t200 * m(3);
t201 = rSges(3,3) + t399;
t404 = t201 * m(3);
t403 = t231 * pkin(2);
t73 = (-t157 * t359 + t317) * t237 + (t185 * t314 + t368) * t231;
t76 = (t157 * t356 + t320) * t237 + t231 * (-t188 * t314 + t371);
t64 = (t238 * t250 - (t251 * t73 + t252 * t76) * t107) * t182;
t402 = t64 * t82;
t74 = (-t157 * t358 + t316) * t239 + (t186 * t313 + t367) * t233;
t77 = (t157 * t355 + t319) * t239 + t233 * (-t189 * t313 + t370);
t65 = (t240 * t250 - (t251 * t74 + t252 * t77) * t108) * t183;
t401 = t65 * t83;
t75 = (-t157 * t357 + t315) * t241 + (t187 * t312 + t366) * t235;
t78 = (t157 * t354 + t318) * t241 + t235 * (-t190 * t312 + t369);
t66 = (t242 * t250 - (t251 * t75 + t252 * t78) * t109) * t184;
t400 = t66 * t84;
t390 = rSges(3,1) * t223;
t383 = t121 * t185;
t382 = t121 * t188;
t381 = t122 * t186;
t380 = t122 * t189;
t379 = t123 * t187;
t378 = t123 * t190;
t212 = sin(t261);
t377 = t124 * t212;
t213 = sin(t262);
t376 = t124 * t213;
t214 = sin(t263);
t375 = t124 * t214;
t158 = 2 * t208;
t142 = sin(t158);
t374 = t134 * t142;
t159 = 2 * t209;
t143 = sin(t159);
t373 = t134 * t143;
t160 = 2 * t211;
t144 = sin(t160);
t372 = t134 * t144;
t145 = cos(t158);
t163 = -rSges(3,1) * t428 + Icges(3,4);
t365 = t163 * t145;
t146 = cos(t159);
t364 = t163 * t146;
t147 = cos(t160);
t363 = t163 * t147;
t164 = rSges(2,1) * t429 - Icges(2,4);
t215 = cos(t261);
t362 = t164 * t215;
t216 = cos(t262);
t361 = t164 * t216;
t217 = cos(t263);
t360 = t164 * t217;
t350 = pkin(3) * t446;
t349 = rSges(2,1) * t417 - Icges(2,5);
t344 = t176 + t223;
t343 = t177 + t223;
t342 = t180 + t223;
t341 = t270 + t271;
t340 = -0.2e1 * t427;
t339 = 0.2e1 * pkin(3);
t338 = -0.2e1 * t413;
t337 = 0.2e1 * t413;
t336 = -0.2e1 * t410;
t335 = -0.2e1 * t407;
t333 = 0.2e1 * t402;
t331 = 0.2e1 * t401;
t329 = 0.2e1 * t400;
t328 = pkin(2) * t406;
t327 = pkin(2) * t405;
t326 = pkin(2) * t404;
t324 = t222 * t428;
t323 = m(3) * t386;
t322 = m(3) * t385;
t321 = m(3) * t384;
t136 = rSges(3,2) * t406 - Icges(3,6);
t139 = rSges(3,1) * t406 - Icges(3,5);
t283 = t243 * t429 - Icges(2,6);
t284 = -t243 * t248 + Icges(2,5);
t67 = (t136 * t222 - t139 * t223 + t284 - t328) * t231 - (t136 * t223 + t139 * t222 + t283) * t237;
t310 = t67 * t386;
t137 = rSges(3,2) * t405 - Icges(3,6);
t140 = rSges(3,1) * t405 - Icges(3,5);
t68 = (t137 * t222 - t140 * t223 + t284 - t327) * t233 - (t137 * t223 + t140 * t222 + t283) * t239;
t309 = t68 * t385;
t138 = rSges(3,2) * t404 - Icges(3,6);
t141 = rSges(3,1) * t404 - Icges(3,5);
t69 = (t138 * t222 - t141 * t223 + t284 - t326) * t235 - (t138 * t223 + t141 * t222 + t283) * t241;
t308 = t69 * t384;
t125 = pkin(3) * t172 + t403;
t307 = t121 * t125 * t79;
t127 = pkin(3) * t173 + t412;
t306 = t122 * t127 * t80;
t126 = pkin(3) * t175 + t411;
t305 = t123 * t126 * t81;
t304 = 0.2e1 * m(2) * t427;
t303 = -0.2e1 * rSges(3,1) * t430;
t302 = t428 * t447;
t162 = m(1) * rSges(1,1) + m(2) * pkin(1);
t301 = -t162 - t430;
t300 = -0.2e1 * pkin(3) * t270 - 0.4e1 * t407;
t293 = t172 * rSges(3,1) + t178 * rSges(3,2);
t292 = t173 * rSges(3,1) + t179 * rSges(3,2);
t291 = t175 * rSges(3,1) + t181 * rSges(3,2);
t117 = g(1) * t188 - g(2) * t185;
t290 = g(3) * t238 + t117 * t232;
t118 = g(1) * t189 - g(2) * t186;
t289 = g(3) * t240 + t118 * t234;
t119 = g(1) * t190 - g(2) * t187;
t288 = g(3) * t242 + t119 * t236;
t110 = -m(3) * t390 - t161 + t324;
t113 = t429 + (rSges(3,1) * t222 + rSges(3,2) * t223) * m(3);
t287 = -t110 * t237 - t113 * t231;
t286 = -t110 * t239 - t113 * t233;
t285 = -t110 * t241 - t113 * t235;
t268 = pkin(5) ^ 2;
t272 = pkin(1) ^ 2;
t279 = -(2 * t268) - 0.2e1 * t272 - t341 + ((-4 * pkin(5) - 2 * pkin(6)) * pkin(6));
t278 = t272 + t264 / 0.2e1 + t266 / 0.2e1 + t224;
t276 = t268 + ((2 * pkin(5) + rSges(2,3)) * rSges(2,3)) + t272 + t265 / 0.2e1 + t267 / 0.2e1;
t153 = pkin(2) * t324;
t275 = Icges(1,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(2,1) / 0.2e1 + Icges(3,1) / 0.2e1 + Icges(2,2) / 0.2e1 + Icges(3,2) / 0.2e1 - t153;
t260 = 2 * pkin(7);
t255 = pkin(1) * g(3);
t230 = xDDP(1);
t229 = xDDP(2);
t228 = xDDP(3);
t166 = pkin(2) * t224 + t410;
t150 = rSges(2,2) * t417 - Icges(2,6);
t149 = t162 * g(3);
t148 = m(1) * rSges(1,2) - t417;
t133 = t148 * g(3);
t131 = t441 + t341;
t116 = g(1) * t187 + g(2) * t190;
t115 = g(1) * t186 + g(2) * t189;
t114 = g(1) * t185 + g(2) * t188;
t94 = -0.2e1 * t153 + (t265 + t267) * m(2) + Icges(2,3) + Icges(3,3) + (t390 * t446 + t264 + t266 + t271) * m(3);
t72 = (t106 * t242 + t100) * t184 * m(3);
t71 = (t105 * t240 + t99) * t183 * m(3);
t70 = (t104 * t238 + t98) * t182 * m(3);
t63 = pkin(1) * t66;
t62 = pkin(1) * t65;
t61 = pkin(1) * t64;
t54 = t147 * t419 + t217 * t420 + t241 * t337 + (t235 * t340 + t276) * m(2) + t163 * t144 - t164 * t214 + (t201 ^ 2 + (-pkin(2) * t174 + t175 * t448) * rSges(3,2) + (t342 * pkin(2) + t181 * t447) * rSges(3,1) + t278) * m(3) + t275;
t53 = t146 * t419 + t216 * t420 + t239 * t337 + (t233 * t340 + t276) * m(2) + t163 * t143 - t164 * t213 + (t200 ^ 2 + (-pkin(2) * t170 + t173 * t448) * rSges(3,2) + (t344 * pkin(2) + t179 * t447) * rSges(3,1) + t278) * m(3) + t275;
t52 = t145 * t419 + t215 * t420 + t237 * t337 + (t231 * t340 + t276) * m(2) + t163 * t142 - t164 * t212 + (t199 ^ 2 + (-pkin(2) * t171 + t172 * t448) * rSges(3,2) + (t343 * pkin(2) + t178 * t447) * rSges(3,1) + t278) * m(3) + t275;
t45 = -t78 * t308 + t94 * t379;
t44 = -t77 * t309 + t94 * t381;
t43 = -t76 * t310 + t94 * t383;
t42 = -t75 * t308 + t94 * t378;
t41 = -t74 * t309 + t94 * t380;
t40 = -t73 * t310 + t94 * t382;
t39 = (t100 * t414 + t242 * t54) * t184;
t38 = (t240 * t53 + t99 * t415) * t183;
t37 = (t238 * t52 + t98 * t416) * t182;
t36 = (t106 * t78 + t60) * t321;
t35 = (t105 * t77 + t58) * t322;
t34 = (t104 * t76 + t56) * t323;
t33 = (t106 * t75 + t59) * t321;
t32 = (t105 * t74 + t57) * t322;
t31 = (t104 * t73 + t55) * t323;
t27 = t69 * t379 - (t60 * t414 + t54 * t78) * t384;
t26 = t68 * t381 - (t58 * t415 + t53 * t77) * t385;
t25 = t67 * t383 - (t56 * t416 + t52 * t76) * t386;
t24 = t69 * t378 - (t59 * t414 + t54 * t75) * t384;
t23 = t68 * t380 - (t57 * t415 + t53 * t74) * t385;
t22 = t67 * t382 - (t55 * t416 + t52 * t73) * t386;
t18 = t63 + t51 / 0.2e1 + t50 / 0.2e1 - t90 / 0.2e1;
t17 = t62 + t49 / 0.2e1 + t48 / 0.2e1 - t89 / 0.2e1;
t16 = t61 + t47 / 0.2e1 + t46 / 0.2e1 - t88 / 0.2e1;
t15 = (-t81 * t131 / (t194 + (t223 * t241 - t351) * pkin(3)) + (t280 * t66 + t443 - t63) * t66) * t184;
t14 = (-t80 * t131 / (t193 + (t223 * t239 - t352) * pkin(3)) + (t281 * t65 + t444 - t62) * t65) * t183;
t13 = (-t79 * t131 / (t192 + (t223 * t237 - t353) * pkin(3)) + (t282 * t64 + t445 - t61) * t64) * t182;
t12 = ((cos((qJ(2,1) - pkin(7))) * t335 + cos((qJ(2,1) + t260)) * t336 + t300 * t181 + t166 * t432 + t442) * t81 * t421 + ((pkin(1) + t346) * t30 - 0.2e1 * t18 * t156 + t443 * t453 + (-t270 * t147 - t271 * t217 + t279 + ((t431 - 2 * qJ(3,1)) * qJ(3,1))) * t66 / 0.2e1 + (-t342 * t66 * t339 + t18 * t432) * t452 + (t126 + (t144 * t270 + t174 * t350 + t214 * t271) * t421) * t84 * t204) * t66) * t184;
t11 = ((cos((qJ(2,2) - pkin(7))) * t335 + cos((qJ(2,2) + t260)) * t336 + t300 * t179 + t166 * t433 + t442) * t80 * t422 + ((pkin(1) + t347) * t29 - 0.2e1 * t17 * t155 + t444 * t453 + (-t146 * t270 - t216 * t271 + t279 + ((t431 - 2 * qJ(3,2)) * qJ(3,2))) * t65 / 0.2e1 + (-t344 * t65 * t339 + t17 * t433) * t452 + (t127 + (t143 * t270 + t170 * t350 + t213 * t271) * t422) * t83 * t203) * t65) * t183;
t10 = ((cos((qJ(2,3) - pkin(7))) * t335 + cos((t260 + qJ(2,3))) * t336 + t300 * t178 + t166 * t434 + t442) * t79 * t423 + ((pkin(1) + t348) * t28 - 0.2e1 * t16 * t154 + t445 * t453 + (-t145 * t270 - t215 * t271 + t279 + ((t431 - 2 * qJ(3,3)) * qJ(3,3))) * t64 / 0.2e1 + (-t343 * t64 * t339 + t16 * t434) * t452 + (t125 + (t142 * t270 + t171 * t350 + t212 * t271) * t423) * t82 * t202) * t64) * t182;
t9 = -t69 * t15 + t94 * t305 + (-t288 * t110 + t113 * t116) * t235 + (t110 * t116 + t288 * t113) * t241 + ((t372 / 0.2e1 - t363 + t375 / 0.2e1 + t360 + (t161 * t235 + t241 * t429) * pkin(1)) * t66 + (t291 * (t63 + 0.2e1 * t51 + 0.2e1 * t50 - 0.2e1 * t90) + (t30 * t435 + t451 * t66) * pkin(2)) * m(3)) * t66;
t8 = -t68 * t14 + t94 * t306 + (-t289 * t110 + t113 * t115) * t233 + (t110 * t115 + t289 * t113) * t239 + ((t373 / 0.2e1 - t364 + t376 / 0.2e1 + t361 + (t161 * t233 + t239 * t429) * pkin(1)) * t65 + (t292 * (t62 + 0.2e1 * t49 + 0.2e1 * t48 - 0.2e1 * t89) + (t29 * t436 + t449 * t65) * pkin(2)) * m(3)) * t65;
t7 = -t67 * t13 + t94 * t307 + (-t290 * t110 + t113 * t114) * t231 + (t110 * t114 + t290 * t113) * t237 + ((t374 / 0.2e1 - t365 + t377 / 0.2e1 + t362 + (t161 * t231 + t237 * t429) * pkin(1)) * t64 + (t293 * (t61 + 0.2e1 * t47 + 0.2e1 * t46 - 0.2e1 * t88) + (t28 * t437 + t450 * t64) * pkin(2)) * m(3)) * t64;
t6 = (-t66 ^ 2 * t201 - g(3) * t236 - t106 * t15 + t119 * t242 - t12 + (t291 + t411) * t329) * m(3);
t5 = (-t65 ^ 2 * t200 - g(3) * t234 - t105 * t14 + t118 * t240 - t11 + (t292 + t412) * t331) * m(3);
t4 = (-t64 ^ 2 * t199 - g(3) * t232 - t104 * t13 + t117 * t238 - t10 + (t293 + t403) * t333) * m(3);
t3 = -t54 * t15 + t69 * t305 - t12 * t414 + t329 * t363 + t66 * t404 * t443 + (-g(3) * t404 + t133 + (-t285 + t301) * t119) * t242 + t236 * ((-t119 * t201 + t255) * m(3) + t149 + t148 * t119 + t285 * g(3)) + (-t372 - t375) * t400 + ((-t141 * t181 + t138 * t175 - (t326 + t349) * t241 + t150 * t235) * t84 + (t303 * t175 - t302 * t181 + t338 * t235 - t304 * t241) * t66) * t84 - 0.2e1 * (t451 * t254 + t360) * t400;
t2 = -t53 * t14 + t68 * t306 - t11 * t415 + t331 * t364 + t65 * t405 * t444 + (-g(3) * t405 + t133 + (-t286 + t301) * t118) * t240 + t234 * ((-t118 * t200 + t255) * m(3) + t149 + t148 * t118 + t286 * g(3)) + (-t373 - t376) * t401 + ((-t140 * t179 + t137 * t173 - (t327 + t349) * t239 + t150 * t233) * t83 + (t303 * t173 - t302 * t179 + t338 * t233 - t304 * t239) * t65) * t83 - 0.2e1 * (t449 * t254 + t361) * t401;
t1 = -t52 * t13 + t67 * t307 - t10 * t416 + t333 * t365 + t64 * t406 * t445 + (-g(3) * t406 + t133 + (-t287 + t301) * t117) * t238 + t232 * ((-t117 * t199 + t255) * m(3) + t149 + t148 * t117 + t287 * g(3)) + (-t374 - t377) * t402 + ((-t139 * t178 + t136 * t172 - (t328 + t349) * t237 + t150 * t231) * t82 + (t303 * t172 - t302 * t178 + t338 * t231 - t304 * t237) * t64) * t82 - 0.2e1 * (t450 * t254 + t362) * t402;
t19 = [t7 * t383 + t8 * t381 + t9 * t379 - m(4) * g(1) + (t45 * t378 + t44 * t380 + t43 * t382) * t229 + (t45 * t379 + t44 * t381 + t43 * t383 + m(4)) * t230 + ((-t100 * t36 + t242 * t27) * t228 - ((t27 * t78 - t36 * t60) * t230 + (t27 * t75 - t36 * t59) * t229 + t78 * t3 + t60 * t6) * t109) * t184 + ((t240 * t26 - t35 * t99) * t228 - ((t26 * t77 - t35 * t58) * t230 + (t26 * t74 - t35 * t57) * t229 + t77 * t2 + t58 * t5) * t108) * t183 + ((t238 * t25 - t34 * t98) * t228 - ((t25 * t76 - t34 * t56) * t230 + (t25 * t73 - t34 * t55) * t229 + t76 * t1 + t56 * t4) * t107) * t182; t7 * t382 + t8 * t380 + t9 * t378 - m(4) * g(2) + (t42 * t379 + t41 * t381 + t40 * t383) * t230 + (t42 * t378 + t41 * t380 + t40 * t382 + m(4)) * t229 + ((-t100 * t33 + t24 * t242) * t228 - ((t24 * t78 - t33 * t60) * t230 + (t24 * t75 - t33 * t59) * t229 + t75 * t3 + t59 * t6) * t109) * t184 + ((t23 * t240 - t32 * t99) * t228 - ((t23 * t77 - t32 * t58) * t230 + (t23 * t74 - t32 * t57) * t229 + t74 * t2 + t57 * t5) * t108) * t183 + ((t22 * t238 - t31 * t98) * t228 - ((t22 * t76 - t31 * t56) * t230 + (t22 * t73 - t31 * t55) * t229 + t73 * t1 + t55 * t4) * t107) * t182; (-g(3) + t228) * m(4) + ((t72 * t228 + t6) * t100 + (t39 * t228 + t3 + (t187 * t230 + t190 * t229) * t69 * t123) * t242 - ((t39 * t78 + t60 * t72) * t230 + (t39 * t75 + t59 * t72) * t229) * t109) * t184 + ((t71 * t228 + t5) * t99 + (t38 * t228 + t2 + (t186 * t230 + t189 * t229) * t68 * t122) * t240 - ((t38 * t77 + t58 * t71) * t230 + (t38 * t74 + t57 * t71) * t229) * t108) * t183 + ((t70 * t228 + t4) * t98 + (t37 * t228 + t1 + (t185 * t230 + t188 * t229) * t67 * t121) * t238 - ((t37 * t76 + t56 * t70) * t230 + (t37 * t73 + t55 * t70) * t229) * t107) * t182;];
tauX  = t19;
