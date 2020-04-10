% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRR1G2P2A0
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
% Datum: 2020-03-09 21:18
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:18:17
% EndTime: 2020-03-09 21:18:19
% DurationCPUTime: 1.14s
% Computational Cost: add. (14082->129), mult. (6774->259), div. (2988->5), fcn. (7908->18), ass. (0->143)
t487 = 0.1e1 / pkin(2);
t474 = pkin(7) + qJ(2,2);
t463 = qJ(3,2) + t474;
t450 = sin(t463);
t453 = cos(t463);
t457 = sin(t474);
t460 = cos(t474);
t547 = 0.1e1 / (-t450 * t460 + t453 * t457);
t553 = t487 * t547;
t473 = pkin(7) + qJ(2,3);
t462 = qJ(3,3) + t473;
t449 = sin(t462);
t452 = cos(t462);
t456 = sin(t473);
t459 = cos(t473);
t548 = 0.1e1 / (-t449 * t459 + t452 * t456);
t552 = t487 * t548;
t551 = t449 * t548;
t550 = t450 * t547;
t475 = pkin(7) + qJ(2,1);
t464 = qJ(3,1) + t475;
t451 = sin(t464);
t454 = cos(t464);
t458 = sin(t475);
t461 = cos(t475);
t546 = 0.1e1 / (-t451 * t461 + t454 * t458);
t549 = t451 * t546;
t476 = legFrame(3,2);
t468 = cos(t476);
t506 = t468 * t551;
t465 = sin(t476);
t512 = t465 * t551;
t481 = xDP(1);
t524 = t481 * t487;
t480 = xDP(2);
t525 = t480 * t487;
t479 = xDP(3);
t526 = t479 * t487;
t538 = t548 * t452;
t397 = -t506 * t524 + t512 * t525 - t526 * t538;
t427 = pkin(2) * t456 + pkin(3) * t449;
t485 = 0.1e1 / pkin(3);
t523 = t485 * t487;
t511 = t548 * t523;
t529 = t468 * t481;
t532 = t465 * t480;
t430 = pkin(2) * t459 + pkin(3) * t452;
t535 = t430 * t479;
t385 = (t535 / 0.2e1 + (t529 / 0.2e1 - t532 / 0.2e1) * t427) * t511 + t397;
t391 = (t535 + (t529 - t532) * t427) * t511;
t388 = t391 + t397;
t484 = pkin(3) ^ 2;
t486 = pkin(2) ^ 2;
t496 = t449 * t456 + t452 * t459;
t490 = t496 * pkin(2);
t522 = 0.2e1 * pkin(2) * pkin(3);
t542 = t388 * t391;
t379 = ((-t496 * t385 * t522 - t388 * t484 - t397 * t486) * t485 * t397 - (pkin(3) + t490) * t542) * t552;
t382 = (pkin(3) * t542 + (t388 * pkin(3) + t397 * t490) * t397) * t552;
t471 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t433 = -rSges(3,1) * t459 - t456 * rSges(3,2);
t434 = rSges(3,1) * t456 - rSges(3,2) * t459;
t493 = (t433 * t452 - t434 * t449) * pkin(2);
t400 = -Icges(3,3) + (-t471 + t493) * m(3);
t439 = -m(3) * t471 - Icges(3,3);
t545 = (t379 * t439 + t382 * t400) * t548;
t477 = legFrame(2,2);
t469 = cos(t477);
t505 = t469 * t550;
t466 = sin(t477);
t510 = t466 * t550;
t537 = t547 * t453;
t398 = -t505 * t524 + t510 * t525 - t526 * t537;
t428 = pkin(2) * t457 + pkin(3) * t450;
t509 = t547 * t523;
t528 = t469 * t481;
t531 = t466 * t480;
t431 = pkin(2) * t460 + pkin(3) * t453;
t534 = t431 * t479;
t386 = (t534 / 0.2e1 + (t528 / 0.2e1 - t531 / 0.2e1) * t428) * t509 + t398;
t392 = (t534 + (t528 - t531) * t428) * t509;
t389 = t392 + t398;
t495 = t450 * t457 + t453 * t460;
t489 = t495 * pkin(2);
t541 = t389 * t392;
t380 = ((-t495 * t386 * t522 - t389 * t484 - t398 * t486) * t485 * t398 - (pkin(3) + t489) * t541) * t553;
t383 = (pkin(3) * t541 + (t389 * pkin(3) + t398 * t489) * t398) * t553;
t435 = -rSges(3,1) * t460 - t457 * rSges(3,2);
t436 = rSges(3,1) * t457 - rSges(3,2) * t460;
t492 = (t435 * t453 - t436 * t450) * pkin(2);
t401 = -Icges(3,3) + (-t471 + t492) * m(3);
t544 = (t380 * t439 + t383 * t401) * t547;
t478 = legFrame(1,2);
t470 = cos(t478);
t504 = t470 * t549;
t467 = sin(t478);
t508 = t467 * t549;
t536 = t546 * t454;
t399 = -t504 * t524 + t508 * t525 - t526 * t536;
t429 = pkin(2) * t458 + pkin(3) * t451;
t507 = t546 * t523;
t527 = t470 * t481;
t530 = t467 * t480;
t432 = pkin(2) * t461 + pkin(3) * t454;
t533 = t432 * t479;
t387 = (t533 / 0.2e1 + (t527 / 0.2e1 - t530 / 0.2e1) * t429) * t507 + t399;
t393 = (t533 + (t527 - t530) * t429) * t507;
t390 = t393 + t399;
t494 = t451 * t458 + t454 * t461;
t488 = t494 * pkin(2);
t539 = t546 * t487;
t540 = t390 * t393;
t381 = ((-t494 * t387 * t522 - t390 * t484 - t399 * t486) * t485 * t399 - (pkin(3) + t488) * t540) * t539;
t384 = (pkin(3) * t540 + (t390 * pkin(3) + t399 * t488) * t399) * t539;
t437 = -rSges(3,1) * t461 - t458 * rSges(3,2);
t438 = rSges(3,1) * t458 - rSges(3,2) * t461;
t491 = (t437 * t454 - t438 * t451) * pkin(2);
t402 = -Icges(3,3) + (-t471 + t491) * m(3);
t543 = (t381 * t439 + t384 * t402) * t546;
t521 = t427 * t545;
t520 = t428 * t544;
t519 = t429 * t543;
t403 = t433 * t449 + t434 * t452;
t518 = t385 * t391 * t403;
t404 = t435 * t450 + t436 * t453;
t517 = t386 * t392 * t404;
t405 = t437 * t451 + t438 * t454;
t516 = t387 * t393 * t405;
t515 = t397 ^ 2 * t403 * t548;
t514 = t398 ^ 2 * t404 * t547;
t513 = t399 ^ 2 * t405 * t546;
t503 = t427 * t515;
t502 = t428 * t514;
t501 = t429 * t513;
t500 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - Icges(2,3) - Icges(3,3);
t499 = -0.2e1 * t548 * t518;
t498 = -0.2e1 * t547 * t517;
t497 = -0.2e1 * t546 * t516;
t455 = t486 + t471;
t375 = t402 * t381 + (t500 + (-t455 + 0.2e1 * t491) * m(3)) * t384;
t374 = t401 * t380 + (t500 + (-t455 + 0.2e1 * t492) * m(3)) * t383;
t373 = t400 * t379 + (t500 + (-t455 + 0.2e1 * t493) * m(3)) * t382;
t1 = [(-t373 * t506 - t374 * t505 - t375 * t504 + (t468 * t521 + t469 * t520 + t470 * t519) * t485) * t487 + (t468 * t449 * t499 + t469 * t450 * t498 + t470 * t451 * t497 + (-t468 * t503 - t469 * t502 - t470 * t501) * t485) * m(3); (t373 * t512 + t374 * t510 + t375 * t508 + (-t465 * t521 - t466 * t520 - t467 * t519) * t485) * t487 + (0.2e1 * t512 * t518 + 0.2e1 * t510 * t517 + 0.2e1 * t508 * t516 + (t465 * t503 + t466 * t502 + t467 * t501) * t485) * m(3); (-t373 * t538 - t374 * t537 - t375 * t536 + (t430 * t545 + t431 * t544 + t432 * t543) * t485) * t487 + (t452 * t499 + t453 * t498 + t454 * t497 + (-t430 * t515 - t431 * t514 - t432 * t513) * t485) * m(3);];
taucX  = t1;
