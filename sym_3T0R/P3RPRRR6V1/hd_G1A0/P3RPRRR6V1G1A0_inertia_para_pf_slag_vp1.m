% Calculate inertia matrix for parallel robot
% P3RPRRR6V1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:32
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRRR6V1G1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:31:52
% EndTime: 2020-08-06 18:31:53
% DurationCPUTime: 1.02s
% Computational Cost: add. (3980->213), mult. (3674->266), div. (204->7), fcn. (1536->68), ass. (0->160)
t503 = 2 * m(3);
t480 = 2 * qJ(3,3);
t459 = sin(t480);
t467 = sin(qJ(3,3));
t557 = 2 * pkin(2);
t398 = 0.1e1 / (t467 * t557 + pkin(3) * t459 + (sin((pkin(7) + qJ(3,3))) + sin((-pkin(7) + qJ(3,3)))) * pkin(1));
t560 = t398 / 0.2e1;
t481 = 2 * qJ(3,2);
t460 = sin(t481);
t468 = sin(qJ(3,2));
t399 = 0.1e1 / (t468 * t557 + pkin(3) * t460 + (sin((pkin(7) + qJ(3,2))) + sin((-pkin(7) + qJ(3,2)))) * pkin(1));
t559 = t399 / 0.2e1;
t482 = 2 * qJ(3,1);
t461 = sin(t482);
t469 = sin(qJ(3,1));
t400 = 0.1e1 / (t469 * t557 + pkin(3) * t461 + (sin((pkin(7) + qJ(3,1))) + sin((-pkin(7) + qJ(3,1)))) * pkin(1));
t558 = t400 / 0.2e1;
t542 = m(3) / 0.2e1;
t556 = Icges(3,2) / 0.2e1;
t457 = (qJ(1,1) + legFrame(1,3));
t447 = pkin(7) + t457;
t431 = cos(t447);
t440 = t482 + t447;
t443 = -2 * qJ(3,1) + t447;
t452 = qJ(3,1) + t457;
t453 = -qJ(3,1) + t457;
t502 = 0.2e1 * pkin(1);
t441 = qJ(3,1) + t447;
t442 = -qJ(3,1) + t447;
t505 = cos(t441) + cos(t442);
t508 = sin(t441) + sin(t442);
t477 = (-pkin(6) - pkin(5));
t546 = -2 * t477;
t382 = t505 * t557 + (cos(t453) + cos(t452)) * t502 + t508 * t546 + (cos(t443) + cos(t440) + 0.2e1 * t431) * pkin(3);
t555 = t382 * t558;
t456 = (qJ(1,2) + legFrame(2,3));
t446 = pkin(7) + t456;
t430 = cos(t446);
t436 = t481 + t446;
t439 = -2 * qJ(3,2) + t446;
t450 = qJ(3,2) + t456;
t451 = -qJ(3,2) + t456;
t437 = qJ(3,2) + t446;
t438 = -qJ(3,2) + t446;
t506 = cos(t437) + cos(t438);
t509 = sin(t437) + sin(t438);
t381 = t506 * t557 + (cos(t451) + cos(t450)) * t502 + t509 * t546 + (cos(t439) + cos(t436) + 0.2e1 * t430) * pkin(3);
t554 = t381 * t559;
t455 = (qJ(1,3) + legFrame(3,3));
t445 = pkin(7) + t455;
t429 = cos(t445);
t432 = t480 + t445;
t435 = -2 * qJ(3,3) + t445;
t448 = qJ(3,3) + t455;
t449 = -qJ(3,3) + t455;
t433 = qJ(3,3) + t445;
t434 = -qJ(3,3) + t445;
t507 = cos(t433) + cos(t434);
t510 = sin(t433) + sin(t434);
t380 = t507 * t557 + (cos(t449) + cos(t448)) * t502 + t510 * t546 + (cos(t435) + cos(t432) + 0.2e1 * t429) * pkin(3);
t553 = t380 * t560;
t428 = sin(t447);
t545 = 2 * t477;
t379 = t508 * t557 + (sin(t453) + sin(t452)) * t502 + t505 * t545 + (sin(t443) + sin(t440) + 0.2e1 * t428) * pkin(3);
t552 = t379 * t558;
t427 = sin(t446);
t378 = t509 * t557 + (sin(t451) + sin(t450)) * t502 + t506 * t545 + (sin(t439) + sin(t436) + 0.2e1 * t427) * pkin(3);
t551 = t378 * t559;
t426 = sin(t445);
t377 = t510 * t557 + (sin(t449) + sin(t448)) * t502 + t507 * t545 + (sin(t435) + sin(t432) + 0.2e1 * t426) * pkin(3);
t550 = t377 * t560;
t472 = cos(qJ(3,1));
t410 = t472 * rSges(3,1) - t469 * rSges(3,2);
t471 = cos(qJ(3,2));
t409 = t471 * rSges(3,1) - t468 * rSges(3,2);
t470 = cos(qJ(3,3));
t408 = t470 * rSges(3,1) - t467 * rSges(3,2);
t549 = -0.2e1 * pkin(1);
t548 = -2 * pkin(2);
t544 = m(3) * rSges(3,2);
t483 = rSges(3,2) ^ 2;
t484 = rSges(3,1) ^ 2;
t543 = (-t483 + t484) * t542 + t556 - Icges(3,1) / 0.2e1;
t541 = pkin(1) * sin(pkin(7));
t458 = cos(pkin(7)) * pkin(1);
t540 = pkin(2) + t458;
t389 = t429 * t546 + sin(t455) * t549 + t426 * t548 - t510 * pkin(3);
t527 = t389 * t398;
t390 = t430 * t546 + sin(t456) * t549 + t427 * t548 - t509 * pkin(3);
t526 = t390 * t399;
t391 = t431 * t546 + sin(t457) * t549 + t428 * t548 - t508 * pkin(3);
t525 = t391 * t400;
t392 = t426 * t545 + cos(t455) * t549 + t429 * t548 - t507 * pkin(3);
t524 = t392 * t398;
t393 = t427 * t545 + cos(t456) * t549 + t430 * t548 - t506 * pkin(3);
t523 = t393 * t399;
t394 = t428 * t545 + cos(t457) * t549 + t431 * t548 - t505 * pkin(3);
t522 = t394 * t400;
t485 = 0.1e1 / pkin(3);
t521 = t398 * t485;
t520 = t399 * t485;
t519 = t400 * t485;
t405 = 0.1e1 / (pkin(3) * t470 + t540);
t518 = t405 * t426;
t517 = t405 * t429;
t406 = 0.1e1 / (pkin(3) * t471 + t540);
t516 = t406 * t427;
t515 = t406 * t430;
t407 = 0.1e1 / (pkin(3) * t472 + t540);
t514 = t407 * t428;
t513 = t407 * t431;
t504 = t483 + t484;
t512 = (m(3) * t504 + Icges(3,3)) * t485;
t511 = m(2) / 0.2e1 + t542;
t490 = m(3) * t408 * t521;
t494 = t398 * t511;
t362 = t377 * t494 + t389 * t490;
t489 = m(3) * t409 * t520;
t493 = t399 * t511;
t363 = t378 * t493 + t390 * t489;
t488 = m(3) * t410 * t519;
t492 = t400 * t511;
t364 = t379 * t492 + t391 * t488;
t365 = t380 * t494 + t392 * t490;
t366 = t381 * t493 + t393 * t489;
t367 = t382 * t492 + t394 * t488;
t500 = t408 * t542;
t499 = t409 * t542;
t498 = t410 * t542;
t473 = rSges(3,3) + pkin(5);
t491 = t473 + t541;
t401 = -t491 * t544 + Icges(3,6);
t402 = m(3) * rSges(3,1) * t491 - Icges(3,5);
t395 = t401 * t470 - t402 * t467;
t497 = t395 * t521;
t396 = t401 * t471 - t402 * t468;
t496 = t396 * t520;
t397 = t401 * t472 - t402 * t469;
t495 = t397 * t519;
t486 = pkin(1) ^ 2;
t487 = Icges(1,3) + Icges(2,3) + ((2 * pkin(2) ^ 2) + 0.2e1 * t473 ^ 2 + 0.2e1 * t486 + t504) * t542 + t473 * t541 * t503 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t556 + Icges(3,1) / 0.2e1 + (t486 + (-0.2e1 * t541 + rSges(2,2)) * rSges(2,2) + (0.2e1 * t458 + rSges(2,1)) * rSges(2,1)) * m(2);
t454 = -rSges(3,1) * t544 + Icges(3,4);
t376 = cos(t482) * t543 + t454 * t461 + (t410 * t458 + (t458 + t410) * pkin(2)) * t503 + t487;
t375 = cos(t481) * t543 + t454 * t460 + (t409 * t458 + (t458 + t409) * pkin(2)) * t503 + t487;
t374 = cos(t480) * t543 + t454 * t459 + (t408 * t458 + (t458 + t408) * pkin(2)) * t503 + t487;
t361 = t376 * t513 + t391 * t495;
t360 = t375 * t515 + t390 * t496;
t359 = t374 * t517 + t389 * t497;
t358 = -t376 * t514 + t394 * t495;
t357 = -t375 * t516 + t393 * t496;
t356 = -t374 * t518 + t392 * t497;
t355 = t397 * t513 + (t379 * t498 + t391 * t512) * t400;
t354 = t396 * t515 + (t378 * t499 + t390 * t512) * t399;
t353 = t395 * t517 + (t377 * t500 + t389 * t512) * t398;
t352 = -t397 * t514 + (t382 * t498 + t394 * t512) * t400;
t351 = -t396 * t516 + (t381 * t499 + t393 * t512) * t399;
t350 = -t395 * t518 + (t380 * t500 + t392 * t512) * t398;
t349 = t367 + t366 + t365;
t348 = t364 + t363 + t362;
t1 = [-t356 * t518 - t357 * t516 - t358 * t514 + m(4) + (t350 * t524 + t351 * t523 + t352 * t522) * t485 + t365 * t553 + t366 * t554 + t367 * t555, t356 * t517 + t357 * t515 + t358 * t513 + (t350 * t527 + t351 * t526 + t352 * t525) * t485 + t365 * t550 + t366 * t551 + t367 * t552, t349; -t359 * t518 - t360 * t516 - t361 * t514 + (t353 * t524 + t354 * t523 + t355 * t522) * t485 + t362 * t553 + t363 * t554 + t364 * t555, t359 * t517 + t360 * t515 + t361 * t513 + m(4) + (t353 * t527 + t354 * t526 + t355 * t525) * t485 + t362 * t550 + t363 * t551 + t364 * t552, t348; t349, t348, 0.3e1 * m(2) + (3 * m(3)) + m(4);];
MX  = t1;
