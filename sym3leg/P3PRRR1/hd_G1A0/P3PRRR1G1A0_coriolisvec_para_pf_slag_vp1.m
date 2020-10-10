% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRR1G1A0
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2020-03-09 21:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRR1G1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR1G1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRR1G1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:14:47
% EndTime: 2020-03-09 21:14:48
% DurationCPUTime: 0.87s
% Computational Cost: add. (8692->121), mult. (4734->234), div. (1380->5), fcn. (5652->18), ass. (0->129)
t412 = pkin(7) + qJ(2,3);
t395 = sin(t412);
t398 = cos(t412);
t372 = -rSges(3,1) * t398 - t395 * rSges(3,2);
t373 = t395 * rSges(3,1) - rSges(3,2) * t398;
t401 = qJ(3,3) + t412;
t388 = sin(t401);
t391 = cos(t401);
t475 = 0.1e1 / (t388 * t398 - t395 * t391);
t477 = (t372 * t388 + t373 * t391) * t475;
t413 = pkin(7) + qJ(2,2);
t396 = sin(t413);
t399 = cos(t413);
t374 = -rSges(3,1) * t399 - t396 * rSges(3,2);
t375 = t396 * rSges(3,1) - rSges(3,2) * t399;
t402 = qJ(3,2) + t413;
t389 = sin(t402);
t392 = cos(t402);
t474 = 0.1e1 / (t389 * t399 - t396 * t392);
t476 = (t374 * t389 + t375 * t392) * t474;
t414 = pkin(7) + qJ(2,1);
t403 = qJ(3,1) + t414;
t390 = sin(t403);
t393 = cos(t403);
t397 = sin(t414);
t400 = cos(t414);
t473 = 0.1e1 / (t390 * t400 - t397 * t393);
t419 = xDP(1);
t425 = 0.1e1 / pkin(2);
t450 = t419 * t425;
t418 = xDP(2);
t451 = t418 * t425;
t417 = legFrame(1,3);
t406 = sin(t417);
t409 = cos(t417);
t371 = -t390 * t406 + t409 * t393;
t454 = t473 * t371;
t370 = t409 * t390 + t406 * t393;
t455 = t473 * t370;
t335 = t450 * t454 + t451 * t455;
t423 = 0.1e1 / pkin(3);
t449 = t423 * t425;
t442 = t473 * t449;
t350 = -pkin(2) * (t397 * t406 - t409 * t400) + t371 * pkin(3);
t460 = t350 * t419;
t347 = pkin(2) * (t397 * t409 + t406 * t400) + t370 * pkin(3);
t463 = t347 * t418;
t426 = (t460 + t463) * t442;
t326 = -t426 + t335;
t472 = pkin(3) * t326;
t415 = legFrame(3,3);
t404 = sin(t415);
t407 = cos(t415);
t367 = -t388 * t404 + t407 * t391;
t458 = t475 * t367;
t366 = t407 * t388 + t404 * t391;
t459 = t475 * t366;
t333 = t450 * t458 + t451 * t459;
t444 = t475 * t449;
t348 = -pkin(2) * (t395 * t404 - t407 * t398) + t367 * pkin(3);
t462 = t348 * t419;
t345 = pkin(2) * (t395 * t407 + t404 * t398) + t366 * pkin(3);
t465 = t345 * t418;
t321 = -(t462 / 0.2e1 + t465 / 0.2e1) * t444 + t333;
t428 = (t462 + t465) * t444;
t324 = -t428 + t333;
t422 = pkin(3) ^ 2;
t424 = pkin(2) ^ 2;
t437 = t388 * t395 + t391 * t398;
t431 = t437 * pkin(2);
t448 = 0.2e1 * pkin(2) * pkin(3);
t453 = t475 * t425;
t468 = t324 * t428;
t315 = ((-t437 * t321 * t448 - t324 * t422 - t424 * t333) * t423 * t333 + (pkin(3) + t431) * t468) * t453;
t318 = (-pkin(3) * t468 + (pkin(3) * t324 + t333 * t431) * t333) * t453;
t410 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t434 = (t372 * t391 - t373 * t388) * pkin(2);
t336 = -Icges(3,3) + (-t410 + t434) * m(3);
t378 = -t410 * m(3) - Icges(3,3);
t471 = (-t378 * t315 - t336 * t318) * t475;
t416 = legFrame(2,3);
t405 = sin(t416);
t408 = cos(t416);
t369 = -t389 * t405 + t408 * t392;
t456 = t474 * t369;
t368 = t408 * t389 + t405 * t392;
t457 = t474 * t368;
t334 = t450 * t456 + t451 * t457;
t443 = t474 * t449;
t349 = -pkin(2) * (t396 * t405 - t408 * t399) + t369 * pkin(3);
t461 = t349 * t419;
t346 = pkin(2) * (t396 * t408 + t405 * t399) + t368 * pkin(3);
t464 = t346 * t418;
t322 = -(t461 / 0.2e1 + t464 / 0.2e1) * t443 + t334;
t427 = (t461 + t464) * t443;
t325 = -t427 + t334;
t436 = t389 * t396 + t392 * t399;
t430 = t436 * pkin(2);
t452 = t474 * t425;
t467 = t325 * t427;
t316 = ((-t436 * t322 * t448 - t325 * t422 - t424 * t334) * t423 * t334 + (pkin(3) + t430) * t467) * t452;
t319 = (-pkin(3) * t467 + (pkin(3) * t325 + t334 * t430) * t334) * t452;
t433 = (t374 * t392 - t375 * t389) * pkin(2);
t337 = -Icges(3,3) + (-t410 + t433) * m(3);
t470 = (-t378 * t316 - t337 * t319) * t474;
t323 = -(t460 / 0.2e1 + t463 / 0.2e1) * t442 + t335;
t435 = t390 * t397 + t393 * t400;
t429 = t435 * pkin(2);
t466 = t426 * t473;
t317 = (-(-t435 * t323 * t448 - t326 * t422 - t424 * t335) * t423 * t473 * t335 - t326 * (pkin(3) + t429) * t466) * t425;
t320 = (-t426 * t472 + (t335 * t429 + t472) * t335) * t425 * t473;
t376 = -rSges(3,1) * t400 - t397 * rSges(3,2);
t377 = t397 * rSges(3,1) - rSges(3,2) * t400;
t432 = (t376 * t393 - t377 * t390) * pkin(2);
t338 = -Icges(3,3) + (-t410 + t432) * m(3);
t469 = (t378 * t317 - t338 * t320) * t473;
t447 = t333 ^ 2 * t477;
t446 = t334 ^ 2 * t476;
t353 = t376 * t390 + t377 * t393;
t445 = t335 ^ 2 * t353 * t473;
t441 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - Icges(2,3) - Icges(3,3);
t440 = -0.2e1 * t321 * t428 * t477;
t439 = -0.2e1 * t322 * t427 * t476;
t438 = -0.2e1 * t323 * t353 * t466;
t394 = t424 + t410;
t311 = t338 * t317 - (t441 + (-t394 + 0.2e1 * t432) * m(3)) * t320;
t310 = -t337 * t316 - (t441 + (-t394 + 0.2e1 * t433) * m(3)) * t319;
t309 = -t336 * t315 - (t441 + (-t394 + 0.2e1 * t434) * m(3)) * t318;
t1 = [(t309 * t458 + t310 * t456 + t311 * t454 + (-t348 * t471 - t349 * t470 - t350 * t469) * t423) * t425 + (t367 * t440 + t369 * t439 + t371 * t438 + (t348 * t447 + t349 * t446 + t350 * t445) * t423) * m(3); (t309 * t459 + t310 * t457 + t311 * t455 + (-t345 * t471 - t346 * t470 - t347 * t469) * t423) * t425 + (t366 * t440 + t368 * t439 + t370 * t438 + (t345 * t447 + t346 * t446 + t347 * t445) * t423) * m(3); 0;];
taucX  = t1;
