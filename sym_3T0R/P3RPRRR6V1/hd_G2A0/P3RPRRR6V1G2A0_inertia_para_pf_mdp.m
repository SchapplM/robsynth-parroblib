% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RPRRR6V1G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRRR6V1G2A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:37
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3RPRRR6V1G2A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:37:32
% EndTime: 2020-08-06 18:37:40
% DurationCPUTime: 7.75s
% Computational Cost: add. (8262->448), mult. (9696->906), div. (1260->20), fcn. (8127->53), ass. (0->374)
t2549 = sin(pkin(7));
t2568 = -pkin(6) - pkin(5);
t2505 = t2549 * t2568 - pkin(1);
t2555 = sin(qJ(1,3));
t2490 = t2505 * t2555;
t2550 = cos(pkin(7));
t2561 = cos(qJ(1,3));
t2704 = t2561 * t2568;
t2707 = t2549 * t2561;
t2813 = t2490 - pkin(2) * t2707 - (pkin(2) * t2555 + t2704) * t2550;
t2557 = sin(qJ(1,2));
t2491 = t2505 * t2557;
t2563 = cos(qJ(1,2));
t2703 = t2563 * t2568;
t2706 = t2549 * t2563;
t2812 = t2491 - pkin(2) * t2706 - (pkin(2) * t2557 + t2703) * t2550;
t2559 = sin(qJ(1,1));
t2492 = t2505 * t2559;
t2565 = cos(qJ(1,1));
t2702 = t2565 * t2568;
t2705 = t2549 * t2565;
t2811 = t2492 - pkin(2) * t2705 - (pkin(2) * t2559 + t2702) * t2550;
t2694 = 2 * MDP(6);
t2810 = MDP(10) / 0.2e1;
t2809 = MDP(11) / 0.2e1;
t2560 = cos(qJ(3,3));
t2515 = t2560 * pkin(3) + pkin(2);
t2524 = t2550 * pkin(1);
t2502 = t2524 + t2515;
t2494 = 0.1e1 / t2502 ^ 2;
t2808 = t2494 / 0.4e1;
t2562 = cos(qJ(3,2));
t2516 = t2562 * pkin(3) + pkin(2);
t2503 = t2524 + t2516;
t2496 = 0.1e1 / t2503 ^ 2;
t2807 = t2496 / 0.4e1;
t2564 = cos(qJ(3,1));
t2517 = t2564 * pkin(3) + pkin(2);
t2504 = t2524 + t2517;
t2498 = 0.1e1 / t2504 ^ 2;
t2806 = t2498 / 0.4e1;
t2497 = 0.1e1 / t2504;
t2536 = qJ(1,1) + pkin(7);
t2523 = sin(t2536);
t2697 = pkin(7) + qJ(3,1);
t2700 = -pkin(7) + qJ(3,1);
t2782 = 0.2e1 * t2568;
t2786 = -0.2e1 * pkin(2);
t2787 = -0.2e1 * pkin(1);
t2468 = t2523 * t2782 + cos(t2536) * t2786 + t2565 * t2787 + (-cos(qJ(1,1) - t2700) - cos(qJ(1,1) + t2697)) * pkin(3);
t2558 = sin(qJ(3,1));
t2783 = 0.2e1 * t2558;
t2474 = pkin(3) * sin(0.2e1 * qJ(3,1)) + pkin(2) * t2783 + (sin(t2697) + sin(t2700)) * pkin(1);
t2765 = t2468 / t2474;
t2658 = t2564 * t2765;
t2611 = t2497 * t2658;
t2805 = t2611 / 0.2e1;
t2659 = t2497 * t2765;
t2612 = t2558 * t2659;
t2804 = t2612 / 0.2e1;
t2495 = 0.1e1 / t2503;
t2535 = qJ(1,2) + pkin(7);
t2522 = sin(t2535);
t2696 = pkin(7) + qJ(3,2);
t2699 = -pkin(7) + qJ(3,2);
t2467 = t2522 * t2782 + cos(t2535) * t2786 + t2563 * t2787 + (-cos(qJ(1,2) - t2699) - cos(qJ(1,2) + t2696)) * pkin(3);
t2556 = sin(qJ(3,2));
t2784 = 0.2e1 * t2556;
t2473 = pkin(3) * sin(0.2e1 * qJ(3,2)) + pkin(2) * t2784 + (sin(t2696) + sin(t2699)) * pkin(1);
t2766 = t2467 / t2473;
t2660 = t2562 * t2766;
t2614 = t2495 * t2660;
t2803 = t2614 / 0.2e1;
t2661 = t2495 * t2766;
t2615 = t2556 * t2661;
t2802 = t2615 / 0.2e1;
t2493 = 0.1e1 / t2502;
t2534 = qJ(1,3) + pkin(7);
t2521 = sin(t2534);
t2695 = pkin(7) + qJ(3,3);
t2698 = -pkin(7) + qJ(3,3);
t2466 = t2521 * t2782 + cos(t2534) * t2786 + t2561 * t2787 + (-cos(qJ(1,3) - t2698) - cos(qJ(1,3) + t2695)) * pkin(3);
t2554 = sin(qJ(3,3));
t2785 = 0.2e1 * t2554;
t2472 = pkin(3) * sin(0.2e1 * qJ(3,3)) + pkin(2) * t2785 + (sin(t2695) + sin(t2698)) * pkin(1);
t2767 = t2466 / t2472;
t2662 = t2560 * t2767;
t2617 = t2493 * t2662;
t2801 = t2617 / 0.2e1;
t2663 = t2493 * t2767;
t2618 = t2554 * t2663;
t2800 = t2618 / 0.2e1;
t2799 = -t2498 * t2523 / 0.2e1;
t2798 = -t2496 * t2522 / 0.2e1;
t2797 = -t2494 * t2521 / 0.2e1;
t2553 = legFrame(1,2);
t2511 = t2553 + t2536;
t2512 = -t2553 + t2536;
t2486 = cos(t2512) + cos(t2511);
t2759 = t2486 * t2497;
t2796 = t2759 / 0.2e1;
t2552 = legFrame(2,2);
t2509 = t2552 + t2535;
t2510 = -t2552 + t2535;
t2485 = cos(t2510) + cos(t2509);
t2760 = t2485 * t2495;
t2795 = t2760 / 0.2e1;
t2551 = legFrame(3,2);
t2507 = t2551 + t2534;
t2508 = -t2551 + t2534;
t2484 = cos(t2508) + cos(t2507);
t2761 = t2484 * t2493;
t2794 = t2761 / 0.2e1;
t2483 = -sin(t2511) + sin(t2512);
t2762 = t2483 * t2497;
t2793 = t2762 / 0.2e1;
t2482 = -sin(t2509) + sin(t2510);
t2763 = t2482 * t2495;
t2792 = t2763 / 0.2e1;
t2481 = -sin(t2507) + sin(t2508);
t2764 = t2481 * t2493;
t2791 = t2764 / 0.2e1;
t2790 = t2481 * t2484;
t2789 = t2482 * t2485;
t2788 = t2483 * t2486;
t2781 = MDP(5) / 0.2e1;
t2780 = MDP(5) / 0.4e1;
t2569 = 0.1e1 / pkin(3);
t2779 = MDP(7) * t2569;
t2778 = MDP(8) * t2569;
t2777 = MDP(9) / pkin(3) ^ 2;
t2719 = t2515 * t2555;
t2460 = (t2704 + t2719) * t2550 - t2490 + t2515 * t2707;
t2525 = sin(t2551);
t2776 = t2460 * t2525;
t2528 = cos(t2551);
t2775 = t2460 * t2528;
t2538 = 0.1e1 / t2554;
t2774 = t2460 * t2538;
t2718 = t2516 * t2557;
t2462 = (t2703 + t2718) * t2550 - t2491 + t2516 * t2706;
t2526 = sin(t2552);
t2773 = t2462 * t2526;
t2529 = cos(t2552);
t2772 = t2462 * t2529;
t2541 = 0.1e1 / t2556;
t2771 = t2462 * t2541;
t2717 = t2517 * t2559;
t2464 = (t2702 + t2717) * t2550 - t2492 + t2517 * t2705;
t2527 = sin(t2553);
t2770 = t2464 * t2527;
t2530 = cos(t2553);
t2769 = t2464 * t2530;
t2544 = 0.1e1 / t2558;
t2768 = t2464 * t2544;
t2758 = t2493 * t2521;
t2757 = t2493 * t2525;
t2756 = t2493 * t2528;
t2755 = t2493 * t2538;
t2754 = t2493 * t2560;
t2752 = t2494 * t2525;
t2751 = t2494 * t2528;
t2537 = t2554 ^ 2;
t2750 = t2494 * t2537;
t2539 = 0.1e1 / t2554 ^ 2;
t2749 = t2494 * t2539;
t2748 = t2494 * t2560;
t2747 = t2495 * t2522;
t2746 = t2495 * t2526;
t2745 = t2495 * t2529;
t2744 = t2495 * t2541;
t2743 = t2495 * t2562;
t2741 = t2496 * t2526;
t2740 = t2496 * t2529;
t2540 = t2556 ^ 2;
t2739 = t2496 * t2540;
t2542 = 0.1e1 / t2556 ^ 2;
t2738 = t2496 * t2542;
t2737 = t2496 * t2562;
t2736 = t2497 * t2523;
t2735 = t2497 * t2527;
t2734 = t2497 * t2530;
t2733 = t2497 * t2544;
t2732 = t2497 * t2564;
t2730 = t2498 * t2527;
t2729 = t2498 * t2530;
t2543 = t2558 ^ 2;
t2728 = t2498 * t2543;
t2545 = 0.1e1 / t2558 ^ 2;
t2727 = t2498 * t2545;
t2726 = t2498 * t2564;
t2513 = pkin(1) * t2549 + pkin(5);
t2725 = t2513 * t2521;
t2724 = t2513 * t2522;
t2723 = t2513 * t2523;
t2722 = t2513 / 0.2e1;
t2721 = t2513 * t2569;
t2514 = t2524 + pkin(2);
t2720 = t2514 / 0.2e1;
t2716 = t2525 * t2554;
t2715 = t2526 * t2556;
t2714 = t2527 * t2558;
t2713 = t2528 * t2554;
t2712 = t2529 * t2556;
t2711 = t2530 * t2558;
t2710 = t2538 * t2560;
t2709 = t2541 * t2562;
t2708 = t2544 * t2564;
t2693 = 0.2e1 * t2569;
t2692 = -0.2e1 * t2720;
t2691 = 0.2e1 * t2720;
t2546 = t2560 ^ 2;
t2690 = pkin(3) * (t2550 * t2555 + t2707) * t2546;
t2547 = t2562 ^ 2;
t2689 = pkin(3) * (t2550 * t2557 + t2706) * t2547;
t2548 = t2564 ^ 2;
t2688 = (t2550 * t2559 + t2705) * t2548 * pkin(3);
t2461 = (t2515 * t2561 - t2555 * t2568) * t2550 - t2505 * t2561 - t2549 * t2719;
t2454 = (-t2461 + t2725) * t2754;
t2687 = MDP(11) * t2454 * t2460;
t2463 = (t2516 * t2563 - t2557 * t2568) * t2550 - t2505 * t2563 - t2549 * t2718;
t2455 = (-t2463 + t2724) * t2743;
t2686 = MDP(11) * t2455 * t2462;
t2465 = (t2517 * t2565 - t2559 * t2568) * t2550 - t2505 * t2565 - t2549 * t2717;
t2456 = (-t2465 + t2723) * t2732;
t2685 = MDP(11) * t2456 * t2464;
t2684 = MDP(11) * t2767;
t2683 = MDP(11) * t2766;
t2682 = MDP(11) * t2765;
t2681 = t2460 ^ 2 * t2749;
t2680 = t2462 ^ 2 * t2738;
t2679 = t2464 ^ 2 * t2727;
t2678 = t2460 * t2752;
t2677 = t2460 * t2751;
t2676 = t2460 * t2721;
t2675 = t2460 * t2710;
t2674 = t2461 * t2538 * t2546;
t2673 = t2462 * t2741;
t2672 = t2462 * t2740;
t2671 = t2462 * t2721;
t2670 = t2462 * t2709;
t2669 = t2463 * t2541 * t2547;
t2668 = t2464 * t2730;
t2667 = t2464 * t2729;
t2666 = t2464 * t2721;
t2665 = t2464 * t2708;
t2664 = t2465 * t2544 * t2548;
t2657 = t2494 * t2790;
t2656 = t2496 * t2789;
t2655 = t2498 * t2788;
t2654 = t2514 * t2758;
t2653 = t2521 * t2750;
t2652 = t2546 * t2749;
t2651 = t2539 * t2748;
t2650 = t2554 * t2748;
t2649 = t2514 * t2747;
t2648 = t2522 * t2739;
t2647 = t2547 * t2738;
t2646 = t2542 * t2737;
t2645 = t2556 * t2737;
t2644 = t2514 * t2736;
t2643 = t2523 * t2728;
t2642 = t2548 * t2727;
t2641 = t2545 * t2726;
t2640 = t2558 * t2726;
t2639 = t2554 * t2722;
t2638 = t2556 * t2722;
t2637 = t2558 * t2722;
t2636 = t2560 * t2722;
t2635 = t2562 * t2722;
t2634 = t2564 * t2722;
t2633 = t2694 / 0.2e1;
t2632 = t2694 / 0.4e1;
t2631 = t2693 / 0.2e1;
t2630 = t2554 * t2692;
t2629 = t2556 * t2692;
t2628 = t2558 * t2692;
t2627 = t2560 * t2691;
t2626 = t2562 * t2691;
t2625 = t2564 * t2691;
t2624 = t2494 * t2675;
t2623 = t2461 * t2651;
t2622 = t2496 * t2670;
t2621 = t2463 * t2646;
t2620 = t2498 * t2665;
t2619 = t2465 * t2641;
t2616 = t2721 * t2767;
t2613 = t2721 * t2766;
t2610 = t2721 * t2765;
t2609 = t2521 * t2650;
t2608 = MDP(10) * t2651;
t2607 = t2522 * t2645;
t2606 = MDP(10) * t2646;
t2605 = t2523 * t2640;
t2604 = MDP(10) * t2641;
t2603 = MDP(7) * t2631;
t2602 = MDP(8) * t2631;
t2601 = t2460 * t2461 * t2652;
t2600 = t2663 * t2774;
t2599 = t2525 * t2624;
t2598 = t2528 * t2624;
t2597 = t2675 * t2721;
t2596 = t2462 * t2463 * t2647;
t2595 = t2661 * t2771;
t2594 = t2526 * t2622;
t2593 = t2529 * t2622;
t2592 = t2670 * t2721;
t2591 = t2464 * t2465 * t2642;
t2590 = t2659 * t2768;
t2589 = t2527 * t2620;
t2588 = t2530 * t2620;
t2587 = t2665 * t2721;
t2445 = -t2525 * t2690 + (pkin(3) * t2713 + t2525 * t2813) * t2560 + t2514 * t2713;
t2446 = t2528 * t2690 + (pkin(3) * t2716 - t2528 * t2813) * t2560 + t2514 * t2716;
t2447 = -t2526 * t2689 + (pkin(3) * t2712 + t2526 * t2812) * t2562 + t2514 * t2712;
t2448 = t2529 * t2689 + (pkin(3) * t2715 - t2529 * t2812) * t2562 + t2514 * t2715;
t2449 = -t2527 * t2688 + (pkin(3) * t2711 + t2527 * t2811) * t2564 + t2514 * t2711;
t2450 = t2530 * t2688 + (pkin(3) * t2714 - t2530 * t2811) * t2564 + t2514 * t2714;
t2571 = pkin(1) ^ 2;
t2574 = t2655 / 0.4e1 + t2656 / 0.4e1 + t2657 / 0.4e1;
t2578 = (-t2483 * t2530 + t2486 * t2527) * t2498 * t2464;
t2579 = (-t2482 * t2529 + t2485 * t2526) * t2496 * t2462;
t2580 = (-t2481 * t2528 + t2484 * t2525) * t2494 * t2460;
t2583 = (t2445 * t2446 * t2749 + t2447 * t2448 * t2738 + t2449 * t2450 * t2727 + t2571 * t2574) * MDP(4) + (-t2525 * t2528 * t2681 - t2526 * t2529 * t2680 - t2527 * t2530 * t2679) * t2777 + (t2640 * t2788 + t2645 * t2789 + t2650 * t2790) * t2632 + (t2537 * t2657 + t2540 * t2656 + t2543 * t2655) * t2780 + MDP(1) * t2574 + ((t2578 * t2708 + t2579 * t2709 + t2580 * t2710) * MDP(8) + (t2578 + t2579 + t2580) * MDP(7)) * t2569 / 0.2e1;
t2573 = t2481 * t2797 + t2482 * t2798 + t2483 * t2799;
t2582 = (t2445 * t2623 + t2447 * t2621 + t2449 * t2619 + t2571 * t2573) * MDP(4) + (t2481 * t2801 + t2482 * t2803 + t2483 * t2805 - t2521 * t2599 - t2522 * t2594 - t2523 * t2589) * t2778 + (t2481 * t2800 + t2482 * t2802 + t2483 * t2804 - t2521 * t2678 - t2522 * t2673 - t2523 * t2668) * t2779 + (t2525 * t2600 + t2526 * t2595 + t2527 * t2590) * t2777 + (-t2481 * t2609 - t2482 * t2607 - t2483 * t2605) * t2633 + (-t2481 * t2653 - t2482 * t2648 - t2483 * t2643) * t2781 + MDP(1) * t2573;
t2572 = t2484 * t2797 + t2485 * t2798 + t2486 * t2799;
t2581 = (t2446 * t2623 + t2448 * t2621 + t2450 * t2619 + t2571 * t2572) * MDP(4) + (t2484 * t2801 + t2485 * t2803 + t2486 * t2805 + t2521 * t2598 + t2522 * t2593 + t2523 * t2588) * t2778 + (t2484 * t2800 + t2485 * t2802 + t2486 * t2804 + t2521 * t2677 + t2522 * t2672 + t2523 * t2667) * t2779 + (-t2528 * t2600 - t2529 * t2595 - t2530 * t2590) * t2777 + (-t2484 * t2609 - t2485 * t2607 - t2486 * t2605) * t2633 + (-t2484 * t2653 - t2485 * t2648 - t2486 * t2643) * t2781 + MDP(1) * t2572;
t2518 = t2521 ^ 2;
t2519 = t2522 ^ 2;
t2520 = t2523 ^ 2;
t2577 = t2494 * t2518 + t2496 * t2519 + t2498 * t2520;
t2475 = t2481 ^ 2;
t2476 = t2482 ^ 2;
t2477 = t2483 ^ 2;
t2576 = t2475 * t2808 + t2476 * t2807 + t2477 * t2806;
t2478 = t2484 ^ 2;
t2479 = t2485 ^ 2;
t2480 = t2486 ^ 2;
t2575 = t2478 * t2808 + t2479 * t2807 + t2480 * t2806;
t2453 = (t2558 * t2723 + t2664) * t2497;
t2452 = (t2556 * t2724 + t2669) * t2495;
t2451 = (t2554 * t2725 + t2674) * t2493;
t2444 = -t2556 * t2613 - 0.2e1 * t2562 * t2649;
t2443 = -t2554 * t2616 - 0.2e1 * t2560 * t2654;
t2442 = -t2558 * t2610 - 0.2e1 * t2564 * t2644;
t2441 = -t2564 * t2610 + t2644 * t2783;
t2440 = -t2562 * t2613 + t2649 * t2784;
t2439 = -t2560 * t2616 + t2654 * t2785;
t2438 = (t2485 * t2626 + t2529 * t2671) * t2495;
t2437 = (t2482 * t2626 - t2526 * t2671) * t2495;
t2436 = (t2484 * t2627 + t2528 * t2676) * t2493;
t2435 = (t2481 * t2627 - t2525 * t2676) * t2493;
t2434 = (t2486 * t2625 + t2530 * t2666) * t2497;
t2433 = (t2483 * t2625 - t2527 * t2666) * t2497;
t2430 = (t2486 * t2628 + t2530 * t2587) * t2497;
t2429 = (t2483 * t2628 - t2527 * t2587) * t2497;
t2428 = (t2485 * t2629 + t2529 * t2592) * t2495;
t2427 = (t2482 * t2629 - t2526 * t2592) * t2495;
t2426 = (t2484 * t2630 + t2528 * t2597) * t2493;
t2425 = (t2481 * t2630 - t2525 * t2597) * t2493;
t2420 = (-t2486 * t2634 - t2450) * t2497;
t2419 = (-t2485 * t2635 - t2448) * t2495;
t2418 = (-t2484 * t2636 - t2446) * t2493;
t2417 = (-t2483 * t2634 - t2449) * t2497;
t2416 = (-t2482 * t2635 - t2447) * t2495;
t2415 = (-t2481 * t2636 - t2445) * t2493;
t2413 = (t2450 * t2708 - t2486 * t2637) * t2497;
t2412 = (t2448 * t2709 - t2485 * t2638) * t2495;
t2411 = (t2446 * t2710 - t2484 * t2639) * t2493;
t2410 = (t2449 * t2708 - t2483 * t2637) * t2497;
t2409 = (t2447 * t2709 - t2482 * t2638) * t2495;
t2408 = (t2445 * t2710 - t2481 * t2639) * t2493;
t1 = [MDP(1) * t2575 + (t2446 ^ 2 * t2749 + t2448 ^ 2 * t2738 + t2450 ^ 2 * t2727 + t2571 * t2575) * MDP(4) + (t2478 * t2750 + t2479 * t2739 + t2480 * t2728) * t2780 + (t2478 * t2650 + t2479 * t2645 + t2480 * t2640) * t2632 + (-t2484 * t2677 - t2485 * t2672 - t2486 * t2667) * t2603 + (-t2484 * t2598 - t2485 * t2593 - t2486 * t2588) * t2602 + (t2528 ^ 2 * t2681 + t2529 ^ 2 * t2680 + t2530 ^ 2 * t2679) * t2777 + (t2434 * t2796 + t2436 * t2794 + t2438 * t2795 + ((-t2413 * t2733 - t2450 * t2641) * t2769 + (-t2412 * t2744 - t2448 * t2646) * t2772 + (-t2411 * t2755 - t2446 * t2651) * t2775) * t2569) * MDP(10) + (t2426 * t2794 + t2428 * t2795 + t2430 * t2796 + ((-t2420 * t2497 + t2450 * t2498) * t2530 * t2768 + (-t2419 * t2495 + t2448 * t2496) * t2529 * t2771 + (-t2418 * t2493 + t2446 * t2494) * t2528 * t2774) * t2569) * MDP(11) + MDP(12); (t2433 * t2759 + t2435 * t2761 + t2437 * t2760) * t2810 + (t2425 * t2761 + t2427 * t2760 + t2429 * t2759) * t2809 + ((t2450 * t2527 * t2604 + (-t2410 * MDP(10) * t2734 + (-t2417 * t2734 - t2450 * t2730) * MDP(11)) * t2544) * t2464 + (t2448 * t2526 * t2606 + (-t2409 * MDP(10) * t2745 + (-t2416 * t2745 - t2448 * t2741) * MDP(11)) * t2541) * t2462 + (t2446 * t2525 * t2608 + (-t2408 * MDP(10) * t2756 + (-t2415 * t2756 - t2446 * t2752) * MDP(11)) * t2538) * t2460) * t2569 + t2583; (t2442 * t2759 + t2443 * t2761 + t2444 * t2760) * t2810 + (t2439 * t2761 + t2440 * t2760 + t2441 * t2759) * t2809 + ((-t2450 * t2682 + ((t2450 * t2658 - t2453 * t2769) * MDP(10) - t2530 * t2685) * t2544) * t2497 + (-t2448 * t2683 + ((t2448 * t2660 - t2452 * t2772) * MDP(10) - t2529 * t2686) * t2541) * t2495 + (-t2446 * t2684 + ((t2446 * t2662 - t2451 * t2775) * MDP(10) - t2528 * t2687) * t2538) * t2493) * t2569 + t2581; (t2434 * t2762 + t2436 * t2764 + t2438 * t2763) * t2810 + (t2426 * t2764 + t2428 * t2763 + t2430 * t2762) * t2809 + ((-t2449 * t2530 * t2604 + (t2413 * MDP(10) * t2735 + (t2420 * t2735 + t2449 * t2729) * MDP(11)) * t2544) * t2464 + (-t2447 * t2529 * t2606 + (t2412 * MDP(10) * t2746 + (t2419 * t2746 + t2447 * t2740) * MDP(11)) * t2541) * t2462 + (-t2445 * t2528 * t2608 + (t2411 * MDP(10) * t2757 + (t2418 * t2757 + t2445 * t2751) * MDP(11)) * t2538) * t2460) * t2569 + t2583; MDP(1) * t2576 + (t2445 ^ 2 * t2749 + t2447 ^ 2 * t2738 + t2449 ^ 2 * t2727 + t2571 * t2576) * MDP(4) + (t2475 * t2750 + t2476 * t2739 + t2477 * t2728) * t2780 + (t2475 * t2650 + t2476 * t2645 + t2477 * t2640) * t2632 + (t2481 * t2678 + t2482 * t2673 + t2483 * t2668) * t2603 + (t2481 * t2599 + t2482 * t2594 + t2483 * t2589) * t2602 + (t2525 ^ 2 * t2681 + t2526 ^ 2 * t2680 + t2527 ^ 2 * t2679) * t2777 + (t2433 * t2793 + t2435 * t2791 + t2437 * t2792 + ((t2410 * t2733 + t2449 * t2641) * t2770 + (t2409 * t2744 + t2447 * t2646) * t2773 + (t2408 * t2755 + t2445 * t2651) * t2776) * t2569) * MDP(10) + (t2425 * t2791 + t2427 * t2792 + t2429 * t2793 + ((t2417 * t2497 - t2449 * t2498) * t2527 * t2768 + (t2416 * t2495 - t2447 * t2496) * t2526 * t2771 + (t2415 * t2493 - t2445 * t2494) * t2525 * t2774) * t2569) * MDP(11) + MDP(12); (t2442 * t2762 + t2443 * t2764 + t2444 * t2763) * t2810 + (t2439 * t2764 + t2440 * t2763 + t2441 * t2762) * t2809 + ((-t2449 * t2682 + ((t2449 * t2658 + t2453 * t2770) * MDP(10) + t2527 * t2685) * t2544) * t2497 + (-t2447 * t2683 + ((t2447 * t2660 + t2452 * t2773) * MDP(10) + t2526 * t2686) * t2541) * t2495 + (-t2445 * t2684 + ((t2445 * t2662 + t2451 * t2776) * MDP(10) + t2525 * t2687) * t2538) * t2493) * t2569 + t2582; (-t2434 * t2736 - t2436 * t2758 - t2438 * t2747) * MDP(10) + (-t2426 * t2758 - t2428 * t2747 - t2430 * t2736) * MDP(11) + ((t2411 * t2767 + t2412 * t2766 + t2413 * t2765 - t2528 * t2601 - t2529 * t2596 - t2530 * t2591) * MDP(10) + (t2418 * t2767 + t2419 * t2766 + t2420 * t2765 + t2461 * t2598 + t2463 * t2593 + t2465 * t2588) * MDP(11)) * t2569 + t2581; (-t2433 * t2736 - t2435 * t2758 - t2437 * t2747) * MDP(10) + (-t2425 * t2758 - t2427 * t2747 - t2429 * t2736) * MDP(11) + ((t2408 * t2767 + t2409 * t2766 + t2410 * t2765 + t2525 * t2601 + t2526 * t2596 + t2527 * t2591) * MDP(10) + (t2415 * t2767 + t2416 * t2766 + t2417 * t2765 - t2461 * t2599 - t2463 * t2594 - t2465 * t2589) * MDP(11)) * t2569 + t2582; t2577 * MDP(1) + (t2461 ^ 2 * t2652 + t2463 ^ 2 * t2647 + t2465 ^ 2 * t2642 + t2577 * t2571) * MDP(4) + (t2518 * t2750 + t2519 * t2739 + t2520 * t2728) * MDP(5) + (t2518 * t2650 + t2519 * t2645 + t2520 * t2640) * t2694 + (t2468 ^ 2 / t2474 ^ 2 + t2467 ^ 2 / t2473 ^ 2 + t2466 ^ 2 / t2472 ^ 2) * t2777 + (-t2442 * t2736 - t2443 * t2758 - t2444 * t2747 + ((t2497 * t2664 + t2453) * t2765 + (t2495 * t2669 + t2452) * t2766 + (t2493 * t2674 + t2451) * t2767) * t2569) * MDP(10) + (-t2439 * t2758 - t2440 * t2747 - t2441 * t2736 + ((-t2465 * t2732 + t2456) * t2765 + (-t2463 * t2743 + t2455) * t2766 + (-t2461 * t2754 + t2454) * t2767) * t2569) * MDP(11) + MDP(12) + ((-t2521 * t2618 - t2522 * t2615 - t2523 * t2612) * MDP(7) + (-t2521 * t2617 - t2522 * t2614 - t2523 * t2611) * MDP(8)) * t2693;];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
