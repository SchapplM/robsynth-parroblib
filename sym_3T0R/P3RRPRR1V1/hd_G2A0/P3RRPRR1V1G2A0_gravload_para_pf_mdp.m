% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR1V1G2A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR1V1G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-04 17:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR1V1G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_mdp: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:08:36
% EndTime: 2022-11-04 17:08:37
% DurationCPUTime: 0.72s
% Computational Cost: add. (609->120), mult. (1023->219), div. (177->7), fcn. (1077->18), ass. (0->113)
t2592 = MDP(3) - MDP(13);
t2585 = MDP(9) + MDP(11);
t2510 = sin(qJ(1,3));
t2516 = cos(qJ(1,3));
t2506 = legFrame(3,2);
t2490 = sin(t2506);
t2493 = cos(t2506);
t2535 = t2493 * g(1) - t2490 * g(2);
t2472 = g(3) * t2510 - t2535 * t2516;
t2503 = pkin(3) + qJ(3,3);
t2496 = 0.1e1 / t2503;
t2591 = t2472 * t2496;
t2512 = sin(qJ(1,2));
t2518 = cos(qJ(1,2));
t2507 = legFrame(2,2);
t2491 = sin(t2507);
t2494 = cos(t2507);
t2534 = t2494 * g(1) - t2491 * g(2);
t2473 = g(3) * t2512 - t2534 * t2518;
t2504 = pkin(3) + qJ(3,2);
t2497 = 0.1e1 / t2504;
t2590 = t2473 * t2497;
t2514 = sin(qJ(1,1));
t2520 = cos(qJ(1,1));
t2508 = legFrame(1,2);
t2492 = sin(t2508);
t2495 = cos(t2508);
t2533 = t2495 * g(1) - t2492 * g(2);
t2474 = g(3) * t2514 - t2533 * t2520;
t2505 = pkin(3) + qJ(3,1);
t2498 = 0.1e1 / t2505;
t2589 = t2474 * t2498;
t2554 = MDP(10) + MDP(12);
t2515 = cos(qJ(2,3));
t2499 = 0.1e1 / t2515;
t2578 = t2490 * t2499;
t2517 = cos(qJ(2,2));
t2500 = 0.1e1 / t2517;
t2577 = t2491 * t2500;
t2519 = cos(qJ(2,1));
t2501 = 0.1e1 / t2519;
t2576 = t2492 * t2501;
t2575 = t2493 * t2499;
t2574 = t2494 * t2500;
t2573 = t2495 * t2501;
t2572 = t2496 * t2499;
t2571 = t2496 * t2516;
t2570 = t2497 * t2500;
t2569 = t2497 * t2518;
t2568 = t2498 * t2501;
t2567 = t2498 * t2520;
t2509 = sin(qJ(2,3));
t2521 = pkin(1) + pkin(2);
t2566 = t2509 * t2521;
t2565 = t2510 * t2515;
t2511 = sin(qJ(2,2));
t2564 = t2511 * t2521;
t2563 = t2512 * t2517;
t2513 = sin(qJ(2,1));
t2562 = t2513 * t2521;
t2561 = t2514 * t2519;
t2560 = t2515 * t2516;
t2559 = t2515 * t2521;
t2558 = t2517 * t2518;
t2557 = t2517 * t2521;
t2556 = t2519 * t2520;
t2555 = t2519 * t2521;
t2478 = -t2490 * t2565 + t2509 * t2493;
t2553 = t2478 * t2572;
t2479 = -t2491 * t2563 + t2511 * t2494;
t2552 = t2479 * t2570;
t2480 = -t2492 * t2561 + t2513 * t2495;
t2551 = t2480 * t2568;
t2481 = t2490 * t2509 + t2493 * t2565;
t2550 = t2481 * t2572;
t2482 = t2491 * t2511 + t2494 * t2563;
t2549 = t2482 * t2570;
t2483 = t2492 * t2513 + t2495 * t2561;
t2548 = t2483 * t2568;
t2463 = g(3) * (pkin(1) * t2565 - t2516 * qJ(3,3)) - t2535 * (pkin(1) * t2560 + t2510 * qJ(3,3));
t2547 = t2463 * t2572;
t2464 = g(3) * (pkin(1) * t2563 - t2518 * qJ(3,2)) - t2534 * (pkin(1) * t2558 + t2512 * qJ(3,2));
t2546 = t2464 * t2570;
t2465 = g(3) * (pkin(1) * t2561 - t2520 * qJ(3,1)) - t2533 * (pkin(1) * t2556 + t2514 * qJ(3,1));
t2545 = t2465 * t2568;
t2544 = t2472 * t2571;
t2543 = t2473 * t2569;
t2542 = t2474 * t2567;
t2541 = t2472 * t2553;
t2540 = t2473 * t2552;
t2539 = t2474 * t2551;
t2538 = t2472 * t2550;
t2537 = t2473 * t2549;
t2536 = t2474 * t2548;
t2532 = -t2503 * t2516 + t2510 * t2559;
t2531 = -t2504 * t2518 + t2512 * t2557;
t2530 = -t2505 * t2520 + t2514 * t2555;
t2529 = g(3) * t2516 + t2535 * t2510;
t2528 = g(3) * t2518 + t2534 * t2512;
t2527 = g(3) * t2520 + t2533 * t2514;
t2487 = t2490 * g(1) + t2493 * g(2);
t2457 = -t2487 * t2515 + t2509 * t2529;
t2488 = t2491 * g(1) + t2494 * g(2);
t2459 = -t2488 * t2517 + t2511 * t2528;
t2489 = t2492 * g(1) + t2495 * g(2);
t2461 = -t2489 * t2519 + t2513 * t2527;
t2502 = 0.1e1 / t2521;
t2525 = (t2457 * t2578 + t2459 * t2577 + t2461 * t2576) * t2502;
t2524 = (t2457 * t2575 + t2459 * t2574 + t2461 * t2573) * t2502;
t2462 = t2489 * t2513 + t2519 * t2527;
t2460 = t2488 * t2511 + t2517 * t2528;
t2458 = t2487 * t2509 + t2515 * t2529;
t1 = [(t2536 + t2537 + t2538) * MDP(2) - g(1) * MDP(15) + t2585 * (t2481 * t2591 + t2482 * t2590 + t2483 * t2589 + t2525) + t2554 * (-t2509 * t2538 - t2511 * t2537 - t2513 * t2536 + (t2458 * t2578 + t2460 * t2577 + t2462 * t2576) * t2502) + t2592 * (t2527 * t2548 + t2528 * t2549 + t2529 * t2550) + (t2483 * t2545 - (t2492 * t2562 + t2530 * t2495) * t2589 + t2482 * t2546 - (t2491 * t2564 + t2531 * t2494) * t2590 + t2481 * t2547 - (t2490 * t2566 + t2532 * t2493) * t2591 + pkin(1) * t2525) * MDP(14); (t2539 + t2540 + t2541) * MDP(2) - g(2) * MDP(15) + t2585 * (t2478 * t2591 + t2479 * t2590 + t2480 * t2589 + t2524) + t2554 * (-t2509 * t2541 - t2511 * t2540 - t2513 * t2539 + (t2458 * t2575 + t2460 * t2574 + t2462 * t2573) * t2502) + t2592 * (t2527 * t2551 + t2528 * t2552 + t2529 * t2553) + (t2480 * t2545 - (-t2530 * t2492 + t2495 * t2562) * t2589 + t2479 * t2546 - (-t2531 * t2491 + t2494 * t2564) * t2590 + t2478 * t2547 - (-t2532 * t2490 + t2493 * t2566) * t2591 + pkin(1) * t2524) * MDP(14); (t2542 + t2543 + t2544) * MDP(2) + ((t2520 * t2465 - (t2514 * t2505 + t2520 * t2555) * t2474) * t2498 + (t2518 * t2464 - (t2512 * t2504 + t2518 * t2557) * t2473) * t2497 + (t2516 * t2463 - (t2510 * t2503 + t2516 * t2559) * t2472) * t2496) * MDP(14) - g(3) * MDP(15) + t2592 * (t2527 * t2567 + t2528 * t2569 + t2529 * t2571) + t2554 * (-t2509 * t2544 - t2511 * t2543 - t2513 * t2542) + t2585 * (t2556 * t2589 + t2558 * t2590 + t2560 * t2591);];
taugX  = t1;
