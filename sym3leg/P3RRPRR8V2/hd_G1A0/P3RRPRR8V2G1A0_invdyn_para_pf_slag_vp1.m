% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR8V2G1A0
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
% Datum: 2020-08-06 21:04
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:02:20
% EndTime: 2020-08-06 21:02:30
% DurationCPUTime: 9.36s
% Computational Cost: add. (31728->602), mult. (39888->947), div. (5556->15), fcn. (27801->62), ass. (0->414)
t218 = cos(pkin(7));
t433 = pkin(3) * t218;
t147 = pkin(2) + t433;
t232 = cos(qJ(2,3));
t490 = t147 * t232;
t234 = cos(qJ(2,2));
t489 = t147 * t234;
t236 = cos(qJ(2,1));
t488 = t147 * t236;
t422 = pkin(5) + qJ(3,1);
t199 = -pkin(6) - t422;
t188 = 0.1e1 / t199;
t244 = xDP(3);
t192 = t236 * pkin(2);
t209 = qJ(2,1) + pkin(7);
t178 = cos(t209);
t435 = pkin(3) * t178;
t125 = t192 + t435;
t116 = 0.1e1 / t125;
t167 = sin(t209);
t230 = sin(qJ(2,1));
t439 = pkin(2) * t230;
t122 = pkin(3) * t167 + t439;
t255 = 0.2e1 * qJ(2,1);
t208 = t255 + pkin(7);
t166 = sin(t208);
t473 = 0.2e1 * pkin(2);
t364 = pkin(3) * t473;
t212 = sin(t255);
t263 = pkin(2) ^ 2;
t365 = t263 * t212;
t150 = 0.2e1 * t209;
t138 = sin(t150);
t261 = pkin(3) ^ 2;
t368 = t261 * t138;
t275 = t166 * t364 + t365 + t368;
t487 = 2 * pkin(1);
t75 = t122 * t487 + t275;
t460 = t75 / 0.2e1;
t321 = t116 * t460;
t246 = xDP(1);
t222 = legFrame(1,3);
t182 = sin(t222);
t185 = cos(t222);
t155 = t192 + pkin(1);
t231 = sin(qJ(1,1));
t237 = cos(qJ(1,1));
t89 = t155 * t231 + t199 * t237;
t92 = t155 * t237 - t199 * t231;
t95 = -t182 * t231 + t185 * t237;
t60 = -t182 * t89 + t185 * t92 + t95 * t435;
t404 = t246 * t60;
t245 = xDP(2);
t98 = t182 * t237 + t185 * t231;
t63 = t182 * t92 + t185 * t89 + t98 * t435;
t407 = t245 * t63;
t33 = (t244 * t321 + t404 + t407) * t188;
t421 = pkin(5) + qJ(3,2);
t198 = -pkin(6) - t421;
t187 = 0.1e1 / t198;
t191 = t234 * pkin(2);
t207 = qJ(2,2) + pkin(7);
t174 = cos(t207);
t436 = pkin(3) * t174;
t124 = t191 + t436;
t113 = 0.1e1 / t124;
t165 = sin(t207);
t228 = sin(qJ(2,2));
t440 = pkin(2) * t228;
t121 = pkin(3) * t165 + t440;
t253 = 0.2e1 * qJ(2,2);
t206 = t253 + pkin(7);
t164 = sin(t206);
t211 = sin(t253);
t366 = t263 * t211;
t149 = 0.2e1 * t207;
t137 = sin(t149);
t369 = t261 * t137;
t276 = t164 * t364 + t366 + t369;
t74 = t121 * t487 + t276;
t461 = t74 / 0.2e1;
t322 = t113 * t461;
t221 = legFrame(2,3);
t181 = sin(t221);
t184 = cos(t221);
t154 = t191 + pkin(1);
t229 = sin(qJ(1,2));
t235 = cos(qJ(1,2));
t88 = t154 * t229 + t198 * t235;
t91 = t154 * t235 - t198 * t229;
t94 = -t181 * t229 + t184 * t235;
t59 = -t181 * t88 + t184 * t91 + t94 * t436;
t405 = t246 * t59;
t97 = t181 * t235 + t184 * t229;
t62 = t181 * t91 + t184 * t88 + t97 * t436;
t408 = t245 * t62;
t32 = (t244 * t322 + t405 + t408) * t187;
t420 = pkin(5) + qJ(3,3);
t197 = -pkin(6) - t420;
t186 = 0.1e1 / t197;
t190 = t232 * pkin(2);
t205 = qJ(2,3) + pkin(7);
t171 = cos(t205);
t437 = pkin(3) * t171;
t123 = t190 + t437;
t110 = 0.1e1 / t123;
t163 = sin(t205);
t226 = sin(qJ(2,3));
t441 = pkin(2) * t226;
t120 = pkin(3) * t163 + t441;
t251 = 0.2e1 * qJ(2,3);
t204 = t251 + pkin(7);
t162 = sin(t204);
t210 = sin(t251);
t367 = t263 * t210;
t148 = 0.2e1 * t205;
t136 = sin(t148);
t370 = t261 * t136;
t277 = t162 * t364 + t367 + t370;
t73 = t120 * t487 + t277;
t462 = t73 / 0.2e1;
t323 = t110 * t462;
t220 = legFrame(3,3);
t180 = sin(t220);
t183 = cos(t220);
t153 = t190 + pkin(1);
t227 = sin(qJ(1,3));
t233 = cos(qJ(1,3));
t87 = t153 * t227 + t197 * t233;
t90 = t153 * t233 - t197 * t227;
t93 = -t180 * t227 + t183 * t233;
t58 = -t180 * t87 + t183 * t90 + t93 * t437;
t406 = t246 * t58;
t96 = t180 * t233 + t183 * t227;
t61 = t180 * t90 + t183 * t87 + t96 * t437;
t409 = t245 * t61;
t31 = (t244 * t323 + t406 + t409) * t186;
t486 = -t244 / 0.2e1;
t111 = 0.1e1 / t123 ^ 2;
t216 = t244 ^ 2;
t485 = t111 * t216;
t114 = 0.1e1 / t124 ^ 2;
t484 = t114 * t216;
t117 = 0.1e1 / t125 ^ 2;
t483 = t117 * t216;
t203 = t218 ^ 2;
t482 = (-t203 + 0.1e1) * t261;
t217 = sin(pkin(7));
t371 = t217 * t230;
t336 = pkin(3) * t371;
t481 = -t336 + t488;
t372 = t217 * t228;
t337 = pkin(3) * t372;
t480 = -t337 + t489;
t373 = t217 * t226;
t338 = pkin(3) * t373;
t479 = -t338 + t490;
t176 = cos(t208);
t478 = rSges(3,1) * t166 + rSges(3,2) * t176;
t173 = cos(t206);
t477 = rSges(3,1) * t164 + rSges(3,2) * t173;
t170 = cos(t204);
t476 = rSges(3,1) * t162 + rSges(3,2) * t170;
t475 = -2 * pkin(1);
t80 = 0.1e1 / t479;
t434 = pkin(3) * t217;
t83 = t147 * t226 + t232 * t434;
t425 = t83 * t80;
t52 = (t244 * t425 + t245 * t96 + t246 * t93) * t186;
t472 = -0.2e1 * t52;
t81 = 0.1e1 / t480;
t84 = t147 * t228 + t234 * t434;
t424 = t84 * t81;
t53 = (t244 * t424 + t245 * t97 + t246 * t94) * t187;
t471 = -0.2e1 * t53;
t82 = 0.1e1 / t481;
t85 = t147 * t230 + t236 * t434;
t423 = t85 * t82;
t54 = (t244 * t423 + t245 * t98 + t246 * t95) * t188;
t470 = -0.2e1 * t54;
t160 = pkin(2) * t433;
t219 = t263 / 0.2e1;
t362 = t160 + t219;
t469 = -0.2e1 * (t203 - 0.1e1 / 0.2e1) * t261 - 0.2e1 * t362;
t468 = -0.4e1 * pkin(1) * (t261 / 0.2e1 + t362);
t262 = pkin(2) * t263;
t438 = pkin(2) * t261;
t467 = -0.2e1 * t262 - 0.4e1 * t438;
t466 = m(3) * pkin(1);
t247 = pkin(2) * m(3);
t465 = -t52 / 0.4e1;
t464 = -t53 / 0.4e1;
t463 = -t54 / 0.4e1;
t243 = m(2) * rSges(2,1);
t459 = m(2) * rSges(2,2);
t458 = m(3) * rSges(3,2);
t77 = -rSges(3,1) * t171 + rSges(3,2) * t163 - t153;
t457 = m(3) * t77;
t78 = -rSges(3,1) * t174 + rSges(3,2) * t165 - t154;
t456 = m(3) * t78;
t79 = -rSges(3,1) * t178 + rSges(3,2) * t167 - t155;
t455 = m(3) * t79;
t454 = rSges(2,2) * pkin(1);
t151 = t243 + t247;
t345 = t217 * t458;
t416 = rSges(3,1) * t218;
t86 = -m(3) * t416 - t151 + t345;
t453 = g(3) * t86;
t257 = rSges(2,2) ^ 2;
t259 = (rSges(2,1) ^ 2);
t119 = m(3) * t263 + ((-t257 + t259) * m(2)) + Icges(2,2) - Icges(2,1);
t452 = t119 / 0.2e1;
t256 = rSges(3,2) ^ 2;
t258 = rSges(3,1) ^ 2;
t128 = m(3) * (-t256 + t258) - Icges(3,1) + Icges(3,2);
t451 = t128 / 0.2e1;
t449 = -0.3e1 / 0.4e1 * t263;
t238 = (rSges(2,3) + pkin(5));
t448 = m(2) * t238;
t447 = m(3) * t186;
t446 = m(3) * t187;
t445 = m(3) * t188;
t194 = (rSges(3,3) + t420);
t444 = m(3) * t194;
t195 = (rSges(3,3) + t421);
t443 = m(3) * t195;
t196 = (rSges(3,3) + t422);
t442 = m(3) * t196;
t432 = pkin(3) * t263;
t49 = pkin(1) * t52;
t19 = -t49 + t31;
t193 = t261 + t263;
t397 = (0.2e1 * t160 + t193) * t216;
t412 = t186 * t80;
t10 = t111 * t186 / (t190 + (t218 * t232 - t373) * pkin(3)) * t397 - (-t19 * t338 + t479 * t31 - (t338 * t472 - t19) * t490 - (-t232 ^ 2 * t469 + t482) * t52) * t52 * t412;
t103 = -g(1) * t180 + g(2) * t183;
t106 = g(1) * t183 + g(2) * t180;
t130 = rSges(3,2) * t444 - Icges(3,6);
t133 = rSges(3,1) * t444 - Icges(3,5);
t102 = t459 + (rSges(3,1) * t217 + rSges(3,2) * t218) * m(3);
t320 = m(1) * rSges(1,1) + m(2) * pkin(1) + t466;
t280 = -t102 * t226 - t86 * t232 + t320;
t403 = t110 * t244;
t335 = t52 * t403;
t311 = -0.2e1 * t335;
t142 = m(1) * rSges(1,2) - t448;
t316 = t142 - t444;
t326 = t110 * t120 * t485;
t327 = t458 * t487;
t328 = -0.2e1 * rSges(3,1) * t466;
t329 = 2 * m(2) * t454;
t139 = cos(t148);
t156 = -rSges(3,1) * t458 + Icges(3,4);
t157 = rSges(2,1) * t459 - Icges(2,4);
t213 = cos(t251);
t146 = pkin(2) * t345;
t273 = Icges(1,3) + ((rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1)) + Icges(2,1) / 0.2e1 + Icges(3,1) / 0.2e1 + Icges(2,2) / 0.2e1 + Icges(3,2) / 0.2e1 - t146;
t264 = pkin(1) ^ 2;
t274 = (pkin(5) ^ 2) + ((2 * pkin(5) + rSges(2,3)) * rSges(2,3)) + t264 + t257 / 0.2e1 + t259 / 0.2e1;
t281 = t264 + t256 / 0.2e1 + t258 / 0.2e1 + t219;
t353 = t232 * t487;
t356 = -2 * t454;
t361 = t170 + t218;
t34 = t139 * t451 + t213 * t452 + t151 * t353 + (t226 * t356 + t274) * m(2) + t156 * t136 - t157 * t210 + ((t194 ^ 2) + (-pkin(2) * t162 + t163 * t475) * rSges(3,2) + (t361 * pkin(2) + t171 * t487) * rSges(3,1) + t281) * m(3) + t273;
t348 = pkin(2) * t444;
t354 = t151 * t475;
t363 = rSges(2,1) * t448 - Icges(2,5);
t386 = t157 * t213;
t389 = t156 * t139;
t393 = (rSges(2,2) * t448 - Icges(2,6)) * t244;
t396 = t128 * t136;
t400 = t119 * t210;
t285 = t238 * t459 - Icges(2,6);
t286 = -t238 * t243 + Icges(2,5);
t55 = (t130 * t217 - t133 * t218 + t286 - t348) * t226 - (t130 * t218 + t133 * t217 + t285) * t232;
t332 = t73 * t403;
t16 = -t49 - (-t406 / 0.2e1 - t409 / 0.2e1 - t332 / 0.4e1) * t186;
t249 = 0.2e1 * pkin(7);
t168 = cos(t249 + qJ(2,3));
t172 = cos(qJ(2,3) - pkin(7));
t250 = 0.3e1 * qJ(2,3);
t260 = pkin(3) * t261;
t305 = 0.3e1 / 0.4e1 * t263;
t306 = 0.3e1 / 0.4e1 * t261;
t313 = -0.2e1 * t260 - 0.4e1 * t432;
t319 = -t110 * t186 / 0.2e1;
t349 = -0.2e1 * t432;
t350 = -0.2e1 * t438;
t358 = -0.6e1 * t263 - 0.3e1 * t261;
t46 = (t197 ^ 2 + t264) * t52;
t67 = t197 * t403;
t7 = (t168 * t350 + t313 * t171 + t172 * t349 + t232 * t467 + t468) * t319 * t485 + (t277 * t111 * t197 * t486 + ((-t260 * cos(0.3e1 * t205) - t262 * cos(t250)) * t465 - (t370 / 0.2e1 + t367 / 0.2e1) * t67 + (-(-cos(t250 + pkin(7)) - t172) * t52 * t305 + (-t31 * t487 + t358 * t465 + t46) * t171) * pkin(3) + ((-t52 * t449 + t46) * t232 - (-cos(t249 + t250) - t168 - 0.2e1 * t232) * t52 * t306 + (-0.2e1 * t16 * t361 - t67 * t162) * pkin(3) - (t361 * pkin(3) + t353) * t31) * pkin(2) + (-t31 / 0.2e1 - t16) * (t139 * t261 + t213 * t263 + t193)) * t110) * t186 * t52;
t1 = -t34 * t10 + t55 * t326 - t7 * t457 + t311 * t389 - t31 * t444 * t472 + (-t280 * t103 + t316 * t106) * t233 + t227 * (t316 * t103 + t280 * t106) - (-t396 - t400) * t335 + (t110 * t393 * t226 + (-t133 * t171 + t130 * t163 - (t348 + t363) * t232) * t403 - (t328 * t163 - t327 * t171 + t354 * t226 - t329 * t232) * t52) * t403 + 0.2e1 * (t247 * t476 + t386) * t335;
t431 = t1 * t186;
t298 = rSges(3,1) * t163 + rSges(3,2) * t171;
t4 = (-t52 ^ 2 * t194 - t77 * t10 + t103 * t233 - t106 * t227 - t7 + (t298 + t441) * t311) * m(3);
t430 = t186 * t4;
t104 = -g(1) * t181 + g(2) * t184;
t107 = g(1) * t184 + g(2) * t181;
t50 = pkin(1) * t53;
t20 = -t50 + t32;
t411 = t187 * t81;
t11 = t114 * t187 / (t191 + (t218 * t234 - t372) * pkin(3)) * t397 - (-t20 * t337 + t480 * t32 - (t337 * t471 - t20) * t489 - (-t234 ^ 2 * t469 + t482) * t53) * t53 * t411;
t131 = rSges(3,2) * t443 - Icges(3,6);
t134 = rSges(3,1) * t443 - Icges(3,5);
t279 = -t102 * t228 - t86 * t234 + t320;
t402 = t113 * t244;
t334 = t53 * t402;
t309 = -0.2e1 * t334;
t315 = t142 - t443;
t325 = t113 * t121 * t484;
t347 = pkin(2) * t443;
t140 = cos(t149);
t214 = cos(t253);
t352 = t234 * t487;
t360 = t173 + t218;
t35 = t140 * t451 + t214 * t452 + t151 * t352 + (t228 * t356 + t274) * m(2) + t156 * t137 - t157 * t211 + ((t195 ^ 2) + (-pkin(2) * t164 + t165 * t475) * rSges(3,2) + (t360 * pkin(2) + t174 * t487) * rSges(3,1) + t281) * m(3) + t273;
t385 = t157 * t214;
t388 = t156 * t140;
t395 = t128 * t137;
t399 = t119 * t211;
t56 = (t131 * t217 - t134 * t218 + t286 - t347) * t228 - (t131 * t218 + t134 * t217 + t285) * t234;
t169 = cos(t249 + qJ(2,2));
t331 = t74 * t402;
t17 = -t50 - (-t405 / 0.2e1 - t408 / 0.2e1 - t331 / 0.4e1) * t187;
t175 = cos(qJ(2,2) - pkin(7));
t252 = 0.3e1 * qJ(2,2);
t318 = -t113 * t187 / 0.2e1;
t47 = (t198 ^ 2 + t264) * t53;
t68 = t198 * t402;
t8 = (t169 * t350 + t313 * t174 + t175 * t349 + t234 * t467 + t468) * t318 * t484 + (t276 * t114 * t198 * t486 + ((-t260 * cos(0.3e1 * t207) - t262 * cos(t252)) * t464 - (t369 / 0.2e1 + t366 / 0.2e1) * t68 + (-(-cos(t252 + pkin(7)) - t175) * t53 * t305 + (-t32 * t487 + t358 * t464 + t47) * t174) * pkin(3) + ((-t53 * t449 + t47) * t234 - (-cos(t249 + t252) - t169 - 0.2e1 * t234) * t53 * t306 + (-t68 * t164 - 0.2e1 * t17 * t360) * pkin(3) - (t360 * pkin(3) + t352) * t32) * pkin(2) + (-t32 / 0.2e1 - t17) * (t140 * t261 + t214 * t263 + t193)) * t113) * t187 * t53;
t2 = -t35 * t11 + t56 * t325 - t8 * t456 + t309 * t388 - t32 * t443 * t471 + (-t279 * t104 + t315 * t107) * t235 + t229 * (t315 * t104 + t279 * t107) - (-t395 - t399) * t334 + (t113 * t393 * t228 + (-t134 * t174 + t131 * t165 - (t347 + t363) * t234) * t402 - (t328 * t165 - t327 * t174 + t354 * t228 - t329 * t234) * t53) * t402 + 0.2e1 * (t477 * t247 + t385) * t334;
t429 = t187 * t2;
t297 = rSges(3,1) * t165 + rSges(3,2) * t174;
t5 = (-t53 ^ 2 * t195 + t104 * t235 - t107 * t229 - t78 * t11 - t8 + (t297 + t440) * t309) * m(3);
t428 = t187 * t5;
t105 = -g(1) * t182 + g(2) * t185;
t108 = g(1) * t185 + g(2) * t182;
t51 = pkin(1) * t54;
t21 = -t51 + t33;
t410 = t188 * t82;
t12 = t117 * t188 / (t192 + (t218 * t236 - t371) * pkin(3)) * t397 - (-t21 * t336 + t481 * t33 - (t336 * t470 - t21) * t488 - (-t236 ^ 2 * t469 + t482) * t54) * t54 * t410;
t132 = rSges(3,2) * t442 - Icges(3,6);
t135 = rSges(3,1) * t442 - Icges(3,5);
t278 = -t102 * t230 - t86 * t236 + t320;
t401 = t116 * t244;
t333 = t54 * t401;
t307 = -0.2e1 * t333;
t314 = t142 - t442;
t324 = t116 * t122 * t483;
t346 = pkin(2) * t442;
t141 = cos(t150);
t215 = cos(t255);
t351 = t236 * t487;
t359 = t176 + t218;
t36 = t141 * t451 + t215 * t452 + t151 * t351 + (t230 * t356 + t274) * m(2) + t156 * t138 - t157 * t212 + ((t196 ^ 2) + (-pkin(2) * t166 + t167 * t475) * rSges(3,2) + (t359 * pkin(2) + t178 * t487) * rSges(3,1) + t281) * m(3) + t273;
t384 = t157 * t215;
t387 = t156 * t141;
t394 = t128 * t138;
t398 = t119 * t212;
t57 = (t132 * t217 - t135 * t218 + t286 - t346) * t230 - (t132 * t218 + t135 * t217 + t285) * t236;
t177 = cos(qJ(2,1) + t249);
t179 = cos(qJ(2,1) - pkin(7));
t330 = t75 * t401;
t18 = -t51 - (-t404 / 0.2e1 - t407 / 0.2e1 - t330 / 0.4e1) * t188;
t254 = 0.3e1 * qJ(2,1);
t317 = -t116 * t188 / 0.2e1;
t48 = (t199 ^ 2 + t264) * t54;
t69 = t199 * t401;
t9 = (t177 * t350 + t313 * t178 + t179 * t349 + t236 * t467 + t468) * t317 * t483 + (t275 * t117 * t199 * t486 + ((-t260 * cos(0.3e1 * t209) - t262 * cos(t254)) * t463 - (t368 / 0.2e1 + t365 / 0.2e1) * t69 + (-(-cos(t254 + pkin(7)) - t179) * t54 * t305 + (-t33 * t487 + t358 * t463 + t48) * t178) * pkin(3) + ((-t54 * t449 + t48) * t236 - (-cos(t249 + t254) - t177 - 0.2e1 * t236) * t54 * t306 + (-t69 * t166 - 0.2e1 * t18 * t359) * pkin(3) - (t359 * pkin(3) + t351) * t33) * pkin(2) + (-t33 / 0.2e1 - t18) * (t141 * t261 + t215 * t263 + t193)) * t116) * t188 * t54;
t3 = -t36 * t12 + t57 * t324 - t9 * t455 + t307 * t387 - t33 * t442 * t470 + (-t278 * t105 + t314 * t108) * t237 + t231 * (t314 * t105 + t278 * t108) - (-t394 - t398) * t333 + (t116 * t393 * t230 + (-t135 * t178 + t132 * t167 - (t346 + t363) * t236) * t401 - (t328 * t167 - t327 * t178 + t354 * t230 - t329 * t236) * t54) * t401 + 0.2e1 * (t478 * t247 + t384) * t333;
t427 = t188 * t3;
t296 = rSges(3,1) * t167 + rSges(3,2) * t178;
t6 = (-t54 ^ 2 * t196 + t105 * t237 - t108 * t231 - t79 * t12 - t9 + (t296 + t439) * t307) * m(3);
t426 = t188 * t6;
t223 = xDDP(3);
t383 = t186 * t223;
t224 = xDDP(2);
t382 = t186 * t224;
t225 = xDDP(1);
t381 = t186 * t225;
t380 = t187 * t223;
t379 = t187 * t224;
t378 = t187 * t225;
t377 = t188 * t223;
t376 = t188 * t224;
t375 = t188 * t225;
t343 = t83 * t412;
t341 = t84 * t411;
t339 = t85 * t410;
t289 = t103 * t227 + t106 * t233;
t288 = t104 * t229 + t107 * t235;
t287 = t105 * t231 + t108 * t237;
t99 = g(3) * t102;
t76 = -0.2e1 * t146 + ((t257 + t259) * m(2)) + Icges(2,3) + Icges(3,3) + (t416 * t473 + t256 + t258 + t263) * m(3);
t45 = (t79 * t423 + t321) * t445;
t44 = (t78 * t424 + t322) * t446;
t43 = (t77 * t425 + t323) * t447;
t42 = (t79 * t98 + t63) * t445;
t41 = (t78 * t97 + t62) * t446;
t40 = (t77 * t96 + t61) * t447;
t39 = (t79 * t95 + t60) * t445;
t38 = (t78 * t94 + t59) * t446;
t37 = (t77 * t93 + t58) * t447;
t30 = (t36 * t98 + t63 * t455) * t188;
t29 = (t35 * t97 + t62 * t456) * t187;
t28 = (t34 * t96 + t61 * t457) * t186;
t27 = (t36 * t95 + t60 * t455) * t188;
t26 = (t35 * t94 + t59 * t456) * t187;
t25 = (t34 * t93 + t58 * t457) * t186;
t15 = -t36 * t339 + (t57 - t75 * t79 * t445 / 0.2e1) * t116;
t14 = -t35 * t341 + (t56 - t74 * t78 * t446 / 0.2e1) * t113;
t13 = -t34 * t343 + (t55 - t73 * t77 * t447 / 0.2e1) * t110;
t22 = [-(-t27 * t95 - t39 * t60) * t375 - (-t27 * t98 - t39 * t63) * t376 - (-t27 * t423 + (-t39 * t460 + t95 * t57) * t116) * t377 - t95 * t427 - t60 * t426 - (-t26 * t94 - t38 * t59) * t378 - (-t26 * t97 - t38 * t62) * t379 - (-t26 * t424 + (-t38 * t461 + t94 * t56) * t113) * t380 - t94 * t429 - t59 * t428 - (-t25 * t93 - t37 * t58) * t381 - (-t25 * t96 - t37 * t61) * t382 - (-t25 * t425 + (-t37 * t462 + t93 * t55) * t110) * t383 - t93 * t431 - t58 * t430 + (-g(1) + t225) * m(4); -(-t30 * t95 - t42 * t60) * t375 - (-t30 * t98 - t42 * t63) * t376 - (-t30 * t423 + (-t42 * t460 + t98 * t57) * t116) * t377 - t98 * t427 - t63 * t426 - (-t29 * t94 - t41 * t59) * t378 - (-t29 * t97 - t41 * t62) * t379 - (-t29 * t424 + (-t41 * t461 + t97 * t56) * t113) * t380 - t97 * t429 - t62 * t428 - (-t28 * t93 - t40 * t58) * t381 - (-t28 * t96 - t40 * t61) * t382 - (-t28 * t425 + (-t40 * t462 + t96 * t55) * t110) * t383 - t96 * t431 - t61 * t430 + (-g(2) + t224) * m(4); -(t15 * t95 - t45 * t60) * t375 - (t15 * t98 - t45 * t63) * t376 - t3 * t339 + t116 * (-t57 * t12 + t76 * t324 + (-t287 * t86 + t99) * t230 + (t287 * t102 + t453) * t236 - (-(t394 / 0.2e1 - t387 + t398 / 0.2e1 + t384 + (t151 * t230 + t236 * t459) * pkin(1)) * t54 + (t296 * (-t51 - (-t330 - 0.2e1 * t404 - 0.2e1 * t407) * t188) + (0.2e1 * t230 * t33 - t478 * t54) * pkin(2)) * m(3)) * t54) + t75 * t6 * t317 - (t14 * t94 - t44 * t59) * t378 - (t14 * t97 - t44 * t62) * t379 - t2 * t341 + t113 * (-t56 * t11 + t76 * t325 + (-t288 * t86 + t99) * t228 + (t288 * t102 + t453) * t234 - (-(t395 / 0.2e1 - t388 + t399 / 0.2e1 + t385 + (t151 * t228 + t234 * t459) * pkin(1)) * t53 + (t297 * (-t50 - (-t331 - 0.2e1 * t405 - 0.2e1 * t408) * t187) + (0.2e1 * t228 * t32 - t477 * t53) * pkin(2)) * m(3)) * t53) + t74 * t5 * t318 - (t13 * t93 - t43 * t58) * t381 - (t13 * t96 - t43 * t61) * t382 - t1 * t343 + t110 * (-t55 * t10 + t76 * t326 + (-t289 * t86 + t99) * t226 + (t289 * t102 + t453) * t232 - (-(t396 / 0.2e1 - t389 + t400 / 0.2e1 + t386 + (t151 * t226 + t232 * t459) * pkin(1)) * t52 + (t298 * (-t49 - (-t332 - 0.2e1 * t406 - 0.2e1 * t409) * t186) + (0.2e1 * t226 * t31 - t476 * t52) * pkin(2)) * m(3)) * t52) + t73 * t4 * t319 - m(4) * g(3) + (-t15 * t339 + (t116 * t76 - (t57 * t423 - t45 * t460) * t188) * t116 - t14 * t341 + (t113 * t76 - (t56 * t424 - t44 * t461) * t187) * t113 - t13 * t343 + (t110 * t76 - (t55 * t425 - t43 * t462) * t186) * t110 + m(4)) * t223;];
tauX  = t22;
