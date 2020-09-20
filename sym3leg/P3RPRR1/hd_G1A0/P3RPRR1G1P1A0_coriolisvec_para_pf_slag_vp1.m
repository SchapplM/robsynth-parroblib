% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
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
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:23
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:23:11
% EndTime: 2020-03-09 21:23:12
% DurationCPUTime: 0.74s
% Computational Cost: add. (4828->145), mult. (3510->257), div. (840->7), fcn. (2484->32), ass. (0->124)
t451 = 2 * m(3);
t386 = pkin(7) + qJ(3,1);
t376 = cos(t386);
t397 = cos(qJ(3,1));
t413 = -pkin(1) * t376 - pkin(2) * t397;
t373 = sin(t386);
t394 = sin(qJ(3,1));
t416 = pkin(1) * t373 + t394 * pkin(2);
t450 = t413 * rSges(3,1) + t416 * rSges(3,2);
t385 = pkin(7) + qJ(3,2);
t375 = cos(t385);
t396 = cos(qJ(3,2));
t414 = -pkin(1) * t375 - pkin(2) * t396;
t372 = sin(t385);
t393 = sin(qJ(3,2));
t417 = pkin(1) * t372 + t393 * pkin(2);
t449 = t414 * rSges(3,1) + t417 * rSges(3,2);
t384 = pkin(7) + qJ(3,3);
t374 = cos(t384);
t395 = cos(qJ(3,3));
t415 = -pkin(1) * t374 - pkin(2) * t395;
t371 = sin(t384);
t392 = sin(qJ(3,3));
t418 = pkin(1) * t371 + t392 * pkin(2);
t448 = t415 * rSges(3,1) + t418 * rSges(3,2);
t380 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t405 = pkin(1) ^ 2;
t383 = pkin(2) ^ 2 + t405;
t388 = cos(pkin(7));
t445 = pkin(1) * t388;
t446 = pkin(1) * sin(pkin(7));
t447 = -(t383 + t380) * m(3) - (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - Icges(1,3) - Icges(2,3) - Icges(3,3) - 0.2e1 * (m(2) * rSges(2,1) + m(3) * pkin(2)) * t445 + 0.2e1 * m(2) * rSges(2,2) * t446 + (-rSges(2,1) ^ 2 - rSges(2,2) ^ 2 - t405) * m(2);
t444 = pkin(2) * t388;
t398 = xDP(2);
t403 = 0.1e1 / pkin(3);
t428 = t398 * t403;
t377 = legFrame(3,3) + qJ(1,3);
t367 = pkin(7) + t377;
t361 = qJ(3,3) + t367;
t355 = sin(t361);
t334 = -pkin(1) * sin(t377) - pkin(2) * sin(t367) - pkin(3) * t355;
t351 = 0.1e1 / t418;
t440 = t334 * t351;
t322 = t428 * t440;
t399 = xDP(1);
t427 = t399 * t403;
t358 = cos(t361);
t337 = -pkin(1) * cos(t377) - pkin(2) * cos(t367) - pkin(3) * t358;
t437 = t337 * t351;
t325 = t427 * t437;
t313 = t325 + t322;
t433 = t351 * t358;
t434 = t351 * t355;
t319 = t398 * t434 + t399 * t433;
t310 = t319 + t313;
t443 = t310 * t313;
t378 = legFrame(2,3) + qJ(1,2);
t368 = pkin(7) + t378;
t362 = qJ(3,2) + t368;
t356 = sin(t362);
t335 = -pkin(1) * sin(t378) - pkin(2) * sin(t368) - pkin(3) * t356;
t352 = 0.1e1 / t417;
t439 = t335 * t352;
t323 = t428 * t439;
t359 = cos(t362);
t338 = -pkin(1) * cos(t378) - pkin(2) * cos(t368) - pkin(3) * t359;
t436 = t338 * t352;
t326 = t427 * t436;
t314 = t326 + t323;
t431 = t352 * t359;
t432 = t352 * t356;
t320 = t398 * t432 + t399 * t431;
t311 = t320 + t314;
t442 = t311 * t314;
t379 = legFrame(1,3) + qJ(1,1);
t369 = pkin(7) + t379;
t363 = qJ(3,1) + t369;
t357 = sin(t363);
t336 = -pkin(1) * sin(t379) - pkin(2) * sin(t369) - pkin(3) * t357;
t353 = 0.1e1 / t416;
t438 = t336 * t353;
t324 = t428 * t438;
t360 = cos(t363);
t339 = -pkin(1) * cos(t379) - pkin(2) * cos(t369) - pkin(3) * t360;
t435 = t339 * t353;
t327 = t427 * t435;
t315 = t327 + t324;
t429 = t353 * t360;
t430 = t353 * t357;
t321 = t398 * t430 + t399 * t429;
t312 = t321 + t315;
t441 = t312 * t315;
t426 = 0.2e1 * pkin(2) * pkin(3);
t425 = 0.2e1 * pkin(1);
t424 = t319 ^ 2 * (pkin(2) * (t392 * rSges(3,1) + rSges(3,2) * t395) + (rSges(3,1) * t371 + rSges(3,2) * t374) * pkin(1)) * t351;
t423 = t320 ^ 2 * (pkin(2) * (t393 * rSges(3,1) + rSges(3,2) * t396) + (rSges(3,1) * t372 + rSges(3,2) * t375) * pkin(1)) * t352;
t422 = t321 ^ 2 * (pkin(2) * (t394 * rSges(3,1) + rSges(3,2) * t397) + (rSges(3,1) * t373 + rSges(3,2) * t376) * pkin(1)) * t353;
t370 = pkin(2) + t445;
t307 = t325 / 0.2e1 + t322 / 0.2e1 + t319;
t349 = rSges(3,1) * t446 + rSges(3,2) * t370;
t350 = rSges(3,1) * t370 - rSges(3,2) * t446;
t421 = -0.2e1 * t307 * t313 * (t349 * t395 + t392 * t350) * t351;
t308 = t326 / 0.2e1 + t323 / 0.2e1 + t320;
t420 = -0.2e1 * t308 * t314 * (t349 * t396 + t393 * t350) * t352;
t309 = t327 / 0.2e1 + t324 / 0.2e1 + t321;
t419 = -0.2e1 * t309 * t315 * (t349 * t397 + t394 * t350) * t353;
t402 = pkin(3) ^ 2;
t364 = -t380 * m(3) - Icges(3,3);
t330 = -Icges(3,3) + (-t380 + t450) * m(3);
t329 = -Icges(3,3) + (-t380 + t449) * m(3);
t328 = -Icges(3,3) + (-t380 + t448) * m(3);
t306 = (-pkin(3) * t441 + (-pkin(3) * t312 + t321 * t413) * t321) * t353;
t305 = (-pkin(3) * t442 + (-pkin(3) * t311 + t320 * t414) * t320) * t352;
t304 = (-pkin(3) * t443 + (-pkin(3) * t310 + t319 * t415) * t319) * t351;
t303 = (t309 * t397 * t426 + t312 * t402 + t321 * t383 + (pkin(3) * t309 * t376 + t321 * t444) * t425) * t403 * t353 * t321 + (t370 * t397 - t394 * t446 + pkin(3)) / (t370 * t394 + t397 * t446) * t441;
t302 = (t308 * t396 * t426 + t311 * t402 + t320 * t383 + (pkin(3) * t308 * t375 + t320 * t444) * t425) * t403 * t352 * t320 + (t370 * t396 - t393 * t446 + pkin(3)) / (t370 * t393 + t396 * t446) * t442;
t301 = (t307 * t395 * t426 + t310 * t402 + t319 * t383 + (pkin(3) * t307 * t374 + t319 * t444) * t425) * t403 * t351 * t319 + (t370 * t395 - t392 * t446 + pkin(3)) / (t370 * t392 + t395 * t446) * t443;
t300 = t364 * t303 + t330 * t306;
t299 = t364 * t302 + t329 * t305;
t298 = t364 * t301 + t328 * t304;
t297 = t330 * t303 + (t450 * t451 + t447) * t306;
t296 = t329 * t302 + (t449 * t451 + t447) * t305;
t295 = t328 * t301 + (t448 * t451 + t447) * t304;
t1 = [t295 * t433 + t296 * t431 + t297 * t429 + (t298 * t437 + t299 * t436 + t300 * t435) * t403 + (t358 * t421 + t359 * t420 + t360 * t419 + (t337 * t424 + t338 * t423 + t339 * t422) * t403) * m(3); t295 * t434 + t296 * t432 + t297 * t430 + (t298 * t440 + t299 * t439 + t300 * t438) * t403 + (t355 * t421 + t356 * t420 + t357 * t419 + (t334 * t424 + t335 * t423 + t336 * t422) * t403) * m(3); 0;];
taucX  = t1;
